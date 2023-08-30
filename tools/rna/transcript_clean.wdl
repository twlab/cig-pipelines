version development

import "wdl/subworkflows/crop_reference_fasta_headers.wdl"
import "wdl/tasks/gzip.wdl"
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
        # If this option is set, TranscriptClean will only output transcripts that are either canonical
        # or that contain annotated noncanonical junctions to the clean SAM and Fasta files at the end
        # of the run.
        Boolean canonical_only = true
        String docker = "encodedcc/long-read-rna-pipeline:v2.1.0"

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
    }
}
