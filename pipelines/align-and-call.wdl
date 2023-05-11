version development

import "../tools/bwa/align.wdl"
import "../tools/bwa/untar_idx.wdl"
import "../tools/gatk/bqsr.wdl"
import "../tools/gatk/haplotype_caller.wdl"
import "../tools/picard/markdup.wdl"
import "../tools/samtools/sort.wdl"
import "../tools/samtools/stat.wdl"

workflow align_and_call {
    meta {
        author: "Eddie Belter"
        version: "0.1"
        description: "Align and Call Variants Pipleine"
    }

    input {
        String name
        Array[File] fastqs  # read fastqs
        File idx            # tarred BWA index
        File known_sites    # vcf 
    }

    call untar_idx.untar_idx as reference { input:
        idx=idx,
    }

    call align.bwa_align { input:
        name=name,
        fastqs=fastqs,
        reference=reference.path,
    }

    call stat.samtools_stat { input:
        bam=bwa_align.bam,
    } 

    call sort.samtools_sort { input:
        bam=bwa_align.bam,
    } 

    call markdup.picard_markdup { input:
        name=name,
        bam=samtools_sort.sorted_bam,
    }

    call bqsr.gatk_bqsr { input:
        bam=picard_markdup.dedup_bam,
        reference=reference.path,
        known_sites=known_sites,
    }
 
    call haplotype_caller.gatk_haplotype_caller as hc { input:
        bam=gatk_bqsr.recal_bam,
        reference=reference.path,
    }

    output {
        File bam = gatk_bqsr.recal_bam
        File vcf = hc.vcf
        File dedup_metrics = picard_markdup.metrics
        File stats = samtools_stat.stats
    }
}

task untar_reference {
    input {
        File idx
		}

    command {
        mkdir ref
        cd ref
        tar -xvf ${idx}
        index_folder=$(ls)
        reference_fasta=$(ls | head -1)
        reference_folder=$(pwd)
        reference_index_path=$reference_folder/$reference_fasta
        cd ..
		}

    output {
        File fasta = glob("ref/*.fasta")[0]
		}

    runtime {
        docker: "ebelter/linux-tk:latest"
        cpu : 1
        memory : "4 GB"
    }
}
