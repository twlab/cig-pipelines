version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/abra2.wdl"
import "wdl/tasks/bed/bedtools.wdl"
import "wdl/tasks/picard/markdup.wdl"
import "wdl/tasks/freebayes.wdl"
import "wdl/tasks/gatk/realigner_target_creator.wdl"
import "wdl/tasks/samtools.wdl"
import "wdl/tasks/vcallers/deepvariant.wdl"

workflow genome_updated_wgs_from_bam {
  meta {
      author: "Eddie Belter"
      version: "2.0"
      description: "Left/Realign Indels then Call Variants with DeepVariant, starting from a sorted BAM"
  }

  input {
    String sample
    File sorted_bam                   # pre-sorted input BAM
    File ref_fasta                    # uncompressed reference FASTA
    File ref_fai                      # FASTA index
    File ref_dict                     # sequence dictionary
    Boolean realign_bam = true
    Int targets_expansion_bases = 160
    # dockers and resources
    String abra2_docker
    Int abra2_cpu
    Int abra2_memory
    String bedtools_docker
    Int bedtools_cpu
    Int bedtools_memory
    String deepvariant_docker
    Int deepvariant_cpu
    Int deepvariant_memory
    String freebayes_docker
    Int freebayes_cpu
    Int freebayes_memory
    String gatk_docker
    Int gatk_cpu
    Int gatk_memory
    String picard_docker
    Int picard_cpu
    Int picard_memory
    String samtools_docker
    Int samtools_cpu
    Int samtools_memory
  }

  RunEnv bedtools_runenv = {
    "docker": bedtools_docker,
    "cpu": bedtools_cpu,
    "memory": bedtools_memory,
    "disks": 20,
  }

  RunEnv abra2_renenv = {
    "docker": abra2_docker,
    "cpu": abra2_cpu,
    "memory": abra2_memory,
    "disks": 20,
  }

  RunEnv dv_runenv = {
    "docker": deepvariant_docker,
    "cpu": deepvariant_cpu,
    "memory": deepvariant_memory,
    "disks": 20,
  }

  RunEnv freebayes_renenv = {
    "docker": freebayes_docker,
    "cpu": freebayes_cpu,
    "memory": freebayes_memory,
    "disks": 20,
  }

  RunEnv gatk_renenv = {
    "docker": gatk_docker,
    "cpu": gatk_cpu,
    "memory": gatk_memory,
    "disks": 20,
  }

  RunEnv picard_runenv = {
    "docker": picard_docker,
    "cpu": picard_cpu,
    "memory": picard_memory,
    "disks": 20,
  }

  RunEnv samtools_runenv = {
    "docker": samtools_docker,
    "cpu": samtools_cpu,
    "memory": samtools_memory,
    "disks": 20,
  }

  if ( realign_bam ) {
    call freebayes.run_left_shift_bam as left_shift { input:
      in_bam_file=sorted_bam,
      in_reference_file=ref_fasta,
      in_reference_index_file=ref_fai,
      runenv=freebayes_renenv,
    }

    call samtools.index as left_shift_index { input:
      bam=left_shift.output_bam_file,
      runenv=samtools_runenv,
    }

    call realigner_target_creator.run_realigner_target_creator as target_creator { input:
      in_bam_file=left_shift.output_bam_file,
      in_bam_index_file=left_shift_index.bai,
      in_reference_file=ref_fasta,
      in_reference_index_file=ref_fai,
      in_reference_dict_file=ref_dict,
      runenv=gatk_renenv,
    }

    call bedtools.run_slop as expand_targets { input:
      bed_file=target_creator.output_target_bed_file,
      reference_fai=ref_fai,
      bases=targets_expansion_bases,
      runenv=bedtools_runenv,
    }

    call abra2.run_realigner as realign { input:
      in_bam_file=left_shift.output_bam_file,
      in_bam_index_file=left_shift_index.bai,
      in_target_bed_file=expand_targets.slopped_bed_file,
      in_reference_file=ref_fasta,
      in_reference_index_file=ref_fai,
      runenv=abra2_renenv,
    }
  }

  call markdup.run_markdup as picard_markdup { input:
    bam=select_first([realign.indel_realigned_bam, sorted_bam]),
    params='--SORTING_COLLECTION_SIZE_RATIO .15',
    runenv=picard_runenv,
  }

  call samtools.index as samtools_index { input:
    bam=picard_markdup.dedup_bam,
    runenv=samtools_runenv,
  }

  call deepvariant.run_deepvariant as dv { input:
    sample=sample,
    bam=picard_markdup.dedup_bam,
    bai=samtools_index.bai,
    ref_fasta=ref_fasta,
    ref_fai=ref_fai,
    ref_dict=ref_dict,
    runenv=dv_runenv,
  }

  call samtools.stat as samtools_stat { input:
    bam=picard_markdup.dedup_bam,
    runenv=samtools_runenv,
  }

  output {
    File bam = picard_markdup.dedup_bam
    File bai = samtools_index.bai
    File bam_stats = samtools_stat.stats
    File vcf = dv.vcf
    File vcf_tvi = dv.vcf_tbi
  }
}
