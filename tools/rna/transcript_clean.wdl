version development

import "wdl/subworkflows/crop_reference_fasta_headers.wdl"
import "wdl/tasks/gzip.wdl"
import "wdl/tasks/talon.wdl"
import "wdl/tasks/transcriptclean.wdl"


workflow long_rna_post_align {

    input {
        # Prefix that gets added into output filenames. Default "my_experiment", can not be empty.
        String experiment_prefix
        # Input bam
        File bam
        # Reference genome. Fasta format, gzipped.
        File reference_genome
        # Annotation file, gtf format, gzipped.
        File annotation
        # Variants file, vcf format, gzipped.
        File? variants
        # Is the data from "pacbio" or "nanopore"
        String input_type
        # Array[String] of prefixes for naming novel discoveries in eventual TALON runs (default = "TALON").
        # If defined, length of this array needs to be equal to number of replicates.
        Array[String] talon_prefixes = []
        # Genome build name, for TALON. This must be in the initial_talon_db
        String genome_build
        # Annotation name, for creating abundance from talon db. This must be in the initial_talon_db
        String annotation_name
        # If this option is set, TranscriptClean will only output transcripts that are either canonical
        # or that contain annotated noncanonical junctions to the clean SAM and Fasta files at the end
        # of the run.
        Boolean canonical_only = true
        String docker = "encodedcc/long-read-rna-pipeline:v2.1.0"
        String singularity = "docker://encodedcc/long-read-rna-pipeline:v2.1.0"

        # Resouces
        Resources small_task_resources = {
           "cpu": 2,
           "memory_gb": 7,
           "disks": "local-disk 50 SSD",
        }
        Resources medium_task_resources = {
           "cpu": 6,
           "memory_gb": 32,
           "disks": "local-disk 100 SSD",
        }
        Resources large_task_resources = {
           "cpu": 16,
           "memory_gb": 60,
           "disks": "local-disk 150 SSD",
        }
        Resources xlarge_task_resources = {
           "cpu": 20,
           "memory_gb": 120,
           "disks": "local-disk 150 SSD",
        }
    }

    RuntimeEnvironment runtime_environment = {
      "docker": docker,
      "singularity": singularity
    }

    call crop_reference_fasta_headers.crop_reference_fasta_headers as clean_reference {
        input:
            reference_fasta=reference_genome,
            resources=small_task_resources,
            runtime_environment=runtime_environment,
    }

    call gzip.gzip as decompressed_gtf {
        input:
            input_file=annotation,
            output_filename="combined_annotation.gtf",
            params={
                "decompress": true,
                "noname": false,
            },
            resources=small_task_resources,
            runtime_environment=runtime_environment,
    }

    call transcriptclean.get_SJs_from_gtf as get_splice_junctions {
        input:
            annotation_gtf=decompressed_gtf.out,
            reference_fasta=clean_reference.decompressed,
            resources=medium_task_resources,
            output_filename="SJs.txt",
            runtime_environment=runtime_environment,
    }

    String talon_prefix = "TALON"

    call init_talon_db { input:
        annotation_gtf=decompressed_gtf.out,
        annotation_name=annotation_name,
        ref_genome_name=genome_build,
        idprefix=talon_prefix,
        output_prefix=experiment_prefix,
        resources=medium_task_resources,
        runtime_environment=runtime_environment,
    }

    call split_bam { input:
        bam=bam,
    }

    scatter (j in range(length(split_bam.sams))) {
        String sam = split_bam.sams[j]
        call transcriptclean { input:
            sam=sam,
            reference_genome=clean_reference.decompressed,
            splice_junctions=get_splice_junctions.splice_junctions,
            variants=variants,
            output_prefix=experiment_prefix,
            canonical_only=canonical_only,
            resources=medium_task_resources,
            runtime_environment=runtime_environment,
        }
    }

    call transcriptclean_merge_logs_and_report { input:
        output_prefix=experiment_prefix,
        transcript_logs=transcriptclean.transcript_log,
        transcript_error_logs=transcriptclean.transcript_error_log,
        runtime_environment=runtime_environment,
        resources=small_task_resources,
    }

    call merge_sams { input:
        prefix=experiment_prefix,
        header=split_bam.header,
        sams=transcriptclean.corrected_sam,
    }

    call talon.talon_label_reads { input:
            input_sam=merge_sams.sam,
            output_bam_filename=experiment_prefix+"_labeled.bam",
            output_sam_filename=experiment_prefix+"_labeled.sam",
            output_tsv_filename=experiment_prefix+"_labeled.tsv",
            reference_genome=clean_reference.decompressed,
            resources=medium_task_resources,
            runtime_environment=runtime_environment,
    }

    call talon { input:
        talon_db=init_talon_db.database,
        sam=talon_label_reads.labeled_sam,
        genome_build=genome_build,
        output_prefix=experiment_prefix,
        platform=input_type,
        resources=large_task_resources,
        runtime_environment=runtime_environment,
    }

    call create_abundance_from_talon_db { input:
        talon_db=talon.talon_db_out,
        annotation_name=annotation_name,
        genome_build=genome_build,
        output_prefix=experiment_prefix,
        idprefix=talon_prefix,
        resources=medium_task_resources,
        runtime_environment=runtime_environment,
    }

    call create_gtf_from_talon_db { input:
        talon_db=talon.talon_db_out,
        annotation_name=annotation_name,
        genome_build=genome_build,
        output_prefix=experiment_prefix,
        resources=medium_task_resources,
        runtime_environment=runtime_environment,
    }
}

task init_talon_db {
    input {
        File annotation_gtf
        Resources resources
        String annotation_name
        String ref_genome_name
        String output_prefix
        String? idprefix
        RuntimeEnvironment runtime_environment
    }

    command {
        talon_initialize_database \
            --f ~{annotation_gtf} \
            --a ~{annotation_name} \
            --g ~{ref_genome_name} \
            ~{"--idprefix " + idprefix} \
            --o ~{output_prefix}

        python3.7 $(which record_init_db_inputs.py) \
            --annotation_name ~{annotation_name} \
            --genome ~{ref_genome_name} \
            --outfile ~{output_prefix}_talon_inputs.json
        }

    output {
        File database = "~{output_prefix}.db"
        File talon_inputs = "~{output_prefix}_talon_inputs.json"
    }

    runtime {
        cpu: resources.cpu
        memory: "~{resources.memory_gb} GB"
        disks: resources.disks
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task split_bam {
    input {
        File bam
        #RuntimeEnvironment runtime_environment
    }

    Int cpu = 4
    Int samtools_cpu = cpu - 1
    Int memory = 24
    Int samtools_mem = memory - 4

    String output_basename = sub(basename(bam), ".bam", "")
    String sorted_bam = output_basename+".sorted.bam"
    command <<<
        set -ex
        echo SORT $(date)
        time samtools sort ~{bam} -O BAM -o ~{sorted_bam} -@~{samtools_cpu} -m~{samtools_mem}G

        echo INDEX $(date)
        samtools index -b -@~{samtools_cpu} ~{sorted_bam}

        echo AUTOSOMAL CHR $(date)
        for i in {1..22}
        do
            c="chr${i}"
            echo Chromosome: ${c}
            time samtools view ~{sorted_bam} ${c} -M -O SAM -o "~{output_basename}.${c}.sam"
        done

        echo MINOR CHR $(date)
        samtools view -H ~{sorted_bam} | grep '@SQ' | awk '{print $2}' | awk -F: '{print $2}' | sort -n > all_chr
        for i in {1..22}; do echo "chr${i}"; done | sort -n > autosomal_chr
        comm -23 all_chr autosomal_chr > minor_chr
        samtools view -H ~{sorted_bam} | grep -f minor_chr | grep '@SQ' | sed 's/@SQ\t//' | while read sq; do sn=$(echo ${sq} | tr " " "\n" | grep SN | awk -F: '{print $2}'); ln=$(echo ${sq} | tr " " "\n" | grep LN | awk -F: '{print $2}'); echo -e "${sn}\t1\t${ln}"; done > minor_chr.bed
        time samtools view ~{sorted_bam} -L minor_chr.bed -M -O SAM -o ~{output_basename}.minor.sam
        
        for f in *.sam
        do
            [[ -s ${f} ]] || rm -f ${f}
        done

        echo HEADER
        samtools view -H ~{sorted_bam} -O SAM -o  ~{output_basename}.header
    >>>

    output {
        String prefix = output_basename
        File header = glob(output_basename+"*.header")[0]
        Array[File] sams = glob(output_basename+".*.sam")
    }

    runtime {
        cpu: cpu #resources.cpu
        memory: "~{memory} GB" #"~{resources.memory_gb} GB"
        docker: "ebelter/samtools:1.15.1"
        #disks: resources.disks
        #docker: runtime_environment.docker
        #singularity: runtime_environment.singularity
    }
}

task merge_sams {
    input {
        String prefix
        File header
        Array[File] sams
        #RuntimeEnvironment runtime_environment
    }

    Int cpu = 2
    Int memory = 12

    String output_sam = prefix+".sam"
    command <<<
       cp ~{header} ~{output_sam}
       for sam in ~{sep(' ', sams)}
       do
            echo Appending ${sam}
            cat ${sam} >> ~{output_sam}
       done
    >>>

    output {
        File sam = output_sam
    }

    runtime {
        cpu: cpu #resources.cpu
        memory: "~{memory} GB" #"~{resources.memory_gb} GB"
        docker: "ebelter/samtools:1.15.1"
        #disks: resources.disks
        #docker: runtime_environment.docker
        #singularity: runtime_environment.singularity
    }
}

task transcriptclean {
    input {
        File sam
        File reference_genome
        File splice_junctions
        File? variants
        Resources resources
        String output_prefix
        Boolean canonical_only
        RuntimeEnvironment runtime_environment
    }

    command { 
        test -f ~{variants} && gzip -cd ~{variants} > variants.vcf
        python3.7 $(which TranscriptClean.py) --sam ~{sam} \
            --genome ~{reference_genome} \
            --spliceJns ~{splice_junctions} \
            ~{if defined(variants) then "--variants variants.vcf" else ""} \
            --maxLenIndel 5 \
            --maxSJOffset 5 \
            -m true \
            -i true \
            --correctSJs true \
            --primaryOnly \
            --outprefix ~{output_prefix} \
            --threads ~{resources.cpu} \
            ~{if canonical_only then "--canonOnly" else ""}

        Rscript $(which generate_report.R) ~{output_prefix}
    }
    #samtools view -S -b ~{output_prefix}_clean.sam > ~{output_prefix}_clean.bam

    output {
        #File corrected_bam = "~{output_prefix}_clean.bam"
        File corrected_sam = "~{output_prefix}_clean.sam"
        File corrected_fasta = "~{output_prefix}_clean.fa"
        File transcript_log = "~{output_prefix}_clean.log"
        File transcript_error_log = "~{output_prefix}_clean.TE.log"
        File report = "~{output_prefix}_report.pdf"
    }

    runtime {
        cpu: "8"#resources.cpu
        memory: "48 GB"#"~{resources.memory_gb} GB"
        disks: resources.disks
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task transcriptclean_merge_logs_and_report {
    input {
        String output_prefix
        Array[File] transcript_logs
        Array[File] transcript_error_logs
        RuntimeEnvironment runtime_environment
        Resources resources
    }

    Int cpu = 4
    Int memory = 60

    String merged_log = output_prefix+"_clean.log"
    String merged_te_log = output_prefix+"_clean.TE.log"
    command <<<
        echo Merge TranscriptClean logs...
        head -n1 ~{transcript_logs[0]} > ~{merged_log}
        for log in ~{sep(' ', transcript_logs)}
        do
            echo Appending ${log}
            tail -n+2 ${log} >> ~{merged_log}
        done
 
        echo Merge TranscriptClean TE logs...
        head -1 ~{transcript_error_logs[0]} > ~{merged_te_log}
        for log in ~{sep(' ', transcript_error_logs)}
        do
            echo Appending ${log}
            tail -n+2 ${log} >> ~{merged_te_log}
        done

        #echo Run TranscriptClean report...
        #Rscript $(which generate_report.R) ~{output_prefix}
    >>>

    output {
        #File report = glob("*.pdf")[0]
        File transcript_log = "~{merged_log}"
        File transcript_error_log = "~{merged_te_log}"
    }

    runtime {
        cpu: cpu #resources.cpu
        memory: "~{memory} GB" #"~{resources.memory_gb} GB"
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task talon {
    input {
        File talon_db
        File sam
        Resources resources
        String genome_build
        String output_prefix
        String platform
        RuntimeEnvironment runtime_environment
    }

    command {
        export TMPDIR=/tmp
        echo ~{output_prefix},~{output_prefix},~{platform},~{sam} > ~{output_prefix}_talon_config.csv
        cp ~{talon_db} ./~{output_prefix}_talon.db
        talon --f ~{output_prefix}_talon_config.csv \
                                    --db ~{output_prefix}_talon.db \
                                    --build ~{genome_build} \
                                    --o ~{output_prefix}
    }

    output {
        File talon_config = "~{output_prefix}_talon_config.csv"
        File talon_log = "~{output_prefix}_QC.log"
        File talon_db_out = "~{output_prefix}_talon.db"
    }

    runtime {
        cpu: "4"#resources.cpu
        #cpu: resources.cpu
        memory: "48 GB"#"~{resources.memory_gb} GB"
        #memory: "~{resources.memory_gb} GB"
        disks: resources.disks
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task create_abundance_from_talon_db {
    input {
        File talon_db
        Resources resources
        String annotation_name
        String genome_build
        String output_prefix
        String idprefix
        RuntimeEnvironment runtime_environment
    }

    command {
        talon_abundance --db=~{talon_db} \
                        -a ~{annotation_name} \
                        --build ~{genome_build} \
                        --o=~{output_prefix}
        python3.7 $(which calculate_number_of_genes_detected.py) --abundance ~{output_prefix}_talon_abundance.tsv \
                                                                 --counts_colname ~{output_prefix} \
                                                                 --idprefix ~{idprefix} \
                                                                 --outfile ~{output_prefix}_number_of_genes_detected.json
    }

    output {
        File talon_abundance = "~{output_prefix}_talon_abundance.tsv"
        File number_of_genes_detected = "~{output_prefix}_number_of_genes_detected.json"
    }

    runtime {
        cpu: resources.cpu
        memory: "~{resources.memory_gb} GB"
        disks: resources.disks
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task create_gtf_from_talon_db {
    input {
        File talon_db
        Resources resources
        String annotation_name
        String genome_build
        String output_prefix
        RuntimeEnvironment runtime_environment
    }

    command {
        talon_create_GTF --db ~{talon_db} \
                         -a ~{annotation_name} \
                         --build ~{genome_build} \
                         --o ~{output_prefix}
        gzip -n ~{output_prefix}_talon.gtf
    }

    output {
        File gtf = "~{output_prefix}_talon.gtf.gz"
    }

    runtime {
        cpu: resources.cpu
        memory: "~{resources.memory_gb} GB"
        disks: resources.disks
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task calculate_spearman {
    input {
        File rep1_abundance
        File rep2_abundance
        Resources resources
        String rep1_idprefix
        String rep2_idprefix
        String output_prefix
        RuntimeEnvironment runtime_environment
    }

    command {
        python3.7 $(which calculate_correlation.py) --rep1_abundance ~{rep1_abundance} \
                                                    --rep2_abundance ~{rep2_abundance} \
                                                    --rep1_idprefix ~{rep1_idprefix} \
                                                    --rep2_idprefix ~{rep2_idprefix} \
                                                    --outfile ~{output_prefix}_spearman.json
    }

    output {
        File spearman = "~{output_prefix}_spearman.json"
    }

    runtime {
        cpu: resources.cpu
        memory: "~{resources.memory_gb} GB"
        disks: resources.disks
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task skipNfirstlines {
    input {
        File input_file
        Resources resources
        String output_fn
        Int lines_to_skip
        RuntimeEnvironment runtime_environment
    }

    command {
        sed 1,~{lines_to_skip}d ~{input_file} > ~{output_fn}
    }

    output {
        File output_file = output_fn
    }

    runtime {
        cpu: resources.cpu
        memory: "~{resources.memory_gb} GB"
        disks: resources.disks
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}
