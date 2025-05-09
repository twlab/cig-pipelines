version development

import "wdl/structs/runtime.wdl"
import "wdl/structs/runenv.wdl"
import "wdl/tools/picard/collect_rna_seq_metrics.wdl"
import "wdl/tools/picard/collect_insert_size_metrics.wdl"
import "wdl/tools/samtools/index.wdl"

workflow rna {

    input {
        # endedness: paired or single
        String endedness
        # fastqs_R1: fastq.gz files for Read1 (only these if single-ended)
        Array[Array[File]] fastqs_R1
        # fastqs_R2: fastq.gz files for Read2 (omit if single-ended) in order
        # corresponding to fastqs_R1
        Array[Array[File]] fastqs_R2 = []
        # bamroot: root name for output bams. For example foo_bar will
        # create foo_bar_genome.bam and foo_bar_anno.bam
        String bamroot
        # strandedness: is the library strand specific (stranded or unstranded)
        String strandedness
        # strandedness_direction (forward, reverse, unstranded)
        String strandedness_direction
        # chrom_sizes: chromosome sizes file
        File chrom_sizes
        # Switch to false to not run kallisto
        Boolean run_kallisto = true
        # index: aligner index archive (tar.gz)
        File align_index
        Int align_ncpus
        Int align_ramGB
        String? align_disk
        # picard collect rna seq metrics
        File annotation    # GTF GZ
        # kallisto
        Int? kallisto_number_of_threads
        Int? kallisto_ramGB
        File? kallisto_index
        Array[Int] kallisto_fragment_length = []
        Array[Float] kallisto_sd_of_fragment_length = []
        String? kallisto_disk
        Int bam_to_signals_ncpus
        Int bam_to_signals_ramGB
        String? bam_to_signals_disk
        # rsem_index: RSEM index archive (tar.gz)
        File rsem_index
        # rnd_seed: random seed used for rsem
        Int rnd_seed = 12345
        Int rsem_ncpus
        Int rsem_ramGB
        String? rsem_disk
        File rna_qc_tr_id_to_gene_type_tsv
        String? mad_qc_disk
        String? rna_qc_disk
        File? reference_genome
        File? reference_transcriptome
        # Usually for ENCODE experiments this would be the usual annotation file + the tRNAs in gtf.gz format
        Array[File] reference_annotations = []
        # These are for internal use, leave undefined
        Int? kallisto_fragment_length_undefined
        Float? kallisto_sd_undefined
        String docker = "mgibio/mgi:bulk-rna"
        String singularity = "docker://encodedcc/rna-seq-pipeline:1.2.4"
    }

    RuntimeEnvironment runtime_environment = {
      "docker": docker,
      "singularity": singularity
    }

    # dummy variable value for the single-ended case
    Array[Array[File]] fastqs_R2_ = if (endedness == "single") then fastqs_R1 else fastqs_R2

    scatter (i in range(length(fastqs_R1))) {
        call align { input:
            endedness=endedness,
            fastqs_R1=fastqs_R1[i],
            fastqs_R2=fastqs_R2_[i],
            index=align_index,
            bamroot="rep"+(i+1)+bamroot,
            ncpus=align_ncpus,
            ramGB=align_ramGB,
            disks=align_disk,
            runtime_environment=runtime_environment,
        }

        call index.index as index_genome_bam { input:
            bam=align.genomebam,
            docker = "mgibio/samtools:1.15.1",
        }

        call samtools_quickcheck as check_anno { input:
            bam=align.annobam,
            ncpus=bam_to_signals_ncpus,
            ramGB=bam_to_signals_ramGB,
            disks=bam_to_signals_disk,
            runtime_environment=runtime_environment,
        }

        call collect_rna_seq_metrics.collect_rna_seq_metrics { input:
            alignments=align.genomebam,
            annotation=annotation,
            docker="mgibio/picard:2.27.4",
        }

        call collect_insert_size_metrics.collect_insert_size_metrics { input:
            alignments=align.genomebam,
            docker="mgibio/picard:2.27.4",
        }


        File genome_alignment = align.genomebam
        File transcriptome_alignment = align.annobam

        call bam_to_signals { input:
            input_bam=genome_alignment,
            chrom_sizes=chrom_sizes,
            strandedness=strandedness,
            bamroot="rep"+(i+1)+bamroot+"_genome",
            ncpus=bam_to_signals_ncpus,
            ramGB=bam_to_signals_ramGB,
            disks=bam_to_signals_disk,
            runtime_environment=runtime_environment,
        }

        call rsem_quant { input:
            rsem_index=rsem_index,
            rnd_seed=rnd_seed,
            anno_bam=transcriptome_alignment,
            endedness=endedness,
            read_strand=strandedness_direction,
            ncpus=rsem_ncpus,
            ramGB=rsem_ramGB,
            disks=rsem_disk,
            runtime_environment=runtime_environment,
        }

        call stringtie { input:
            sample=bamroot,
            bam=align.genomebam,
            annotation=annotation,
            docker="ctomlins/stringtie2:latest",
            cpus="8",
            memory="20",
        }
    }

    if (run_kallisto) {
      scatter (i in range(length(fastqs_R1))) {
          Float? kallisto_sd = if (length(kallisto_sd_of_fragment_length) > 0) then kallisto_sd_of_fragment_length[i] else kallisto_sd_undefined
          Int? kallisto_fl = if (length(kallisto_fragment_length) > 0) then kallisto_fragment_length[i] else kallisto_fragment_length_undefined
          call kallisto { input:
              fastqs_R1=fastqs_R1[i],
              fastqs_R2=fastqs_R2_[i],
              endedness=endedness,
              strandedness_direction=strandedness_direction,
              kallisto_index=select_first([kallisto_index]),
              number_of_threads=select_first([kallisto_number_of_threads]),
              ramGB=select_first([kallisto_ramGB]),
              fragment_length=kallisto_fl,
              sd_of_fragment_length=kallisto_sd,
              disks=kallisto_disk,
              out_prefix="rep"+(i+1)+bamroot,
              runtime_environment=runtime_environment,
          }
      }
    }

    # if there are exactly two replicates, calculate the madQC metrics and draw a plot
    if (length(fastqs_R1) == 2) {
        call mad_qc { input:
            quants1=rsem_quant.genes_results[0],
            quants2=rsem_quant.genes_results[1],
            disks=mad_qc_disk,
            runtime_environment=runtime_environment,
        }
    }

    scatter (i in range(length(align.annobam))) {
        call rna_qc { input:
            input_bam=align.anno_bam,
            tr_id_to_gene_type_tsv=rna_qc_tr_id_to_gene_type_tsv,
            output_filename="rep"+(i+1)+bamroot+"_qc.json",
            disks=rna_qc_disk,
            runtime_environment=runtime_environment,
        }
    }
}

task align {
    input {
        Array[File] fastqs_R1
        Array[File] fastqs_R2
        String endedness
        File index
        String bamroot
        Int ncpus
        Int ramGB
        String? disks
        RuntimeEnvironment runtime_environment
    }

    command {
        python3 $(which align.py) \
            --fastqs_R1 ~{sep=' ' fastqs_R1} \
            --fastqs_R2 ~{sep=' ' fastqs_R2} \
            --endedness ~{endedness} \
            --index ~{index} \
            ~{"--bamroot " + bamroot} \
            ~{"--ncpus " + ncpus} \
            ~{"--ramGB " + ramGB}
    }

    output {
        File genomebam = "~{bamroot}_genome.bam"
        File annobam = "~{bamroot}_anno.bam"
        File genome_flagstat = "~{bamroot}_genome_flagstat.txt"
        File anno_flagstat = "~{bamroot}_anno_flagstat.txt"
        File log = "~{bamroot}_Log.final.out"
        File genome_flagstat_json = "~{bamroot}_genome_flagstat.json"
        File anno_flagstat_json = "~{bamroot}_anno_flagstat.json"
        File log_json = "~{bamroot}_Log.final.json"
        File python_log = "align.log"
        File splice_juntions = "SJ.out.tab"
    }

    runtime {
      cpu: ncpus
      memory: "~{ramGB} GB"
      disks : select_first([disks,"local-disk 100 SSD"])
      docker: runtime_environment.docker
      singularity: runtime_environment.singularity
    }
}

task samtools_quickcheck {
    input {
        File bam
        Int ncpus
        Int ramGB
        String? disks
        RuntimeEnvironment runtime_environment
    }

    command {
        samtools quickcheck ~{bam}
    }

    runtime {
      cpu: ncpus
      memory: "~{ramGB} GB"
      disks : select_first([disks,"local-disk 100 SSD"])
      docker: runtime_environment.docker
      singularity: runtime_environment.singularity
    }
}

task  bam_to_signals {
    input {
        File? null
        File input_bam
        File chrom_sizes
        String strandedness
        String bamroot
        Int ncpus
        Int ramGB
        String? disks
        RuntimeEnvironment runtime_environment
    }

    command {
        python3 $(which bam_to_signals.py) \
            --bamfile ~{input_bam} \
            --chrom_sizes ~{chrom_sizes} \
            --strandedness ~{strandedness} \
            --bamroot ~{bamroot}
    }

    output {
        File? unique_unstranded = if (strandedness == "unstranded") then glob("*_genome_uniq.bw")[0] else null
        File? all_unstranded = if (strandedness == "unstranded") then glob("*_genome_all.bw")[0] else null
        File? unique_plus = if (strandedness == "stranded") then glob("*_genome_plusUniq.bw")[0] else null
        File? unique_minus = if (strandedness == "stranded") then glob("*_genome_minusUniq.bw")[0] else null
        File? all_plus = if (strandedness == "stranded") then glob("*_genome_plusAll.bw")[0] else null
        File? all_minus = if (strandedness == "stranded") then glob("*_genome_minusAll.bw")[0] else null
        File python_log = "bam_to_signals.log"
    }

    runtime {
        cpu: ncpus
        memory: "~{ramGB} GB"
        disks : select_first([disks,"local-disk 100 SSD"])
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task rsem_quant {
    input {
        File rsem_index
        File anno_bam
        String endedness
        String read_strand
        Int rnd_seed
        Int ncpus
        Int ramGB
        String? disks
        RuntimeEnvironment runtime_environment
    }

    command {
        python3 $(which rsem_quant.py) \
            --rsem_index ~{rsem_index} \
            --anno_bam ~{anno_bam} \
            --endedness ~{endedness} \
            --read_strand ~{read_strand} \
            --rnd_seed ~{rnd_seed} \
            --ncpus ~{ncpus} \
            --ramGB ~{ramGB}
    }

    output {
        File genes_results = glob("*.genes.results")[0]
        File isoforms_results = glob("*.isoforms.results")[0]
        File python_log = "rsem_quant.log"
        File number_of_genes = glob("*_number_of_genes_detected.json")[0]
    }

    runtime {
        cpu: ncpus
        memory: "~{ramGB} GB"
        disks : select_first([disks,"local-disk 100 SSD"])
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task stringtie {
    input {
        String sample
        File bam
        File annotation
        String docker
        String cpus
        String memory
    }

    String annotation_gunzip = sub(basename("~{annotation}"), ".gz$", "")
    String denovo_transcipts_fn = "~{sample}.stringtie2.denovo.gtf"
    String abundance_estimate_fn = "~{sample}.stringtie2.abundance_estimate.gtf"
    command <<<
        # gunzip anno
        gunzip -c ~{annotation} > ~{annotation_gunzip}
        # denovo
        stringtie -o ~{denovo_transcipts_fn} -p ~{cpus} -m 100 ~{bam}
        # quant
        stringtie -o ~{abundance_estimate_fn} -p ~{cpus} -m 100 -e -b ~{sample} -G ~{annotation_gunzip} ~{bam}
    >>>

    output {
        File denovo_transcipts = glob("~{denovo_transcipts_fn}")[0]
        File abundance_estimate = glob("~{abundance_estimate_fn}")[0]
    }

    runtime {
        docker: docker
        cpu: cpus
        memory: "~{memory} GB"
    }
}

task kallisto {
    input {
        Array[File] fastqs_R1
        Array[File] fastqs_R2
        File kallisto_index
        String endedness
        String strandedness_direction
        Int number_of_threads
        Int ramGB
        String out_prefix
        Int? fragment_length
        Float? sd_of_fragment_length
        String? disks
        RuntimeEnvironment runtime_environment
    }

    command {
        python3 $(which kallisto_quant.py) \
            --fastqs_R1 ~{sep=' ' fastqs_R1} \
            --fastqs_R2 ~{sep=' ' fastqs_R2} \
            --number_of_threads ~{number_of_threads} \
            --strandedness ~{strandedness_direction} \
            --path_to_index ~{kallisto_index} \
            --endedness ~{endedness} \
            ~{"--fragment_length " + fragment_length} \
            ~{"--sd_of_fragment_length " + sd_of_fragment_length} \
            ~{"--out_prefix " + out_prefix}
    }

    output {
        File quants = "kallisto_out/~{out_prefix}_abundance.tsv"
        File python_log = "kallisto_quant.log"
    }

    runtime {
        cpu: number_of_threads
        memory: "~{ramGB} GB"
        disks: select_first([disks, "local-disk 100 SSD"])
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task mad_qc {
    input {
        File quants1
        File quants2
        String? disks
        RuntimeEnvironment runtime_environment
    }

    command {
        python3 $(which mad_qc.py) \
            --quants1 ~{quants1} \
            --quants2 ~{quants2} \
            --MAD_R_path $(which MAD.R)
    }

    output {
        File madQCplot = glob("*_mad_plot.png")[0]
        File madQCmetrics = glob("*_mad_qc_metrics.json")[0]
        File python_log = "mad_qc.log"
    }

    runtime {
        cpu: 2
        memory: "3400 MB"
        disks: select_first([disks,"local-disk 100 SSD"])
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task rna_qc {
    input {
        File input_bam
        File tr_id_to_gene_type_tsv
        String output_filename
        String? disks
        RuntimeEnvironment runtime_environment
    }

    command {
        python3 $(which rna_qc.py) \
            --input_bam ~{input_bam} \
            --tr_id_to_gene_type_tsv ~{tr_id_to_gene_type_tsv} \
            --output_filename ~{output_filename}
    }

    output {
        File rnaQC = output_filename
        File python_log = "rna_qc.log"
    }

    runtime {
        cpu: 2
        memory: "1024 MB"
        disks: select_first([disks, "local-disk 100 SSD"])
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}
