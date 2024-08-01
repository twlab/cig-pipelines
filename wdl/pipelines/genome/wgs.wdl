version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/abra2.wdl"
import "wdl/tasks/bed/bedtools.wdl"
import "wdl/tasks/bwa/align.wdl"
import "wdl/tasks/bwa/idx.wdl"
import "wdl/tasks/freebayes.wdl"
import "wdl/tasks/gatk/realigner_target_creator.wdl"
import "wdl/tasks/samtools.wdl"
import "wdl/tasks/vcallers/deepvariant.wdl"

workflow genome_wgs_with_realign {
  meta {
      author: "Eddie Belter"
      version: "1.2"
      description: "Align with BWA, Left/Realign Indels then Call Variants with DeepVariant"
  }

  input {
    String sample
    Array[File] fastqs
    File idx # tarred BWA index
    Boolean realign_bam = true
    Int targets_expansion_bases = 160
    # dockers
    String abra2_docker = "mgibio/abra2:v2.23-focal"
    String bwa_docker = "ebelter/bwa:0.7.17"
    String bedtools_docker = "biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1"
    String deepvariant_docker = "google/deepvariant:1.5.0"
    String freebayes_docker = "mgibio/freebayes:1.3.6-focal"
    String gatk_docker = "broadinstitute/gatk3:3.5-0"
    String samtools_docker = "mgibio/samtools:1.15.1-buster"
  }

  # RunEnvs in order of usage
  RunEnv idx_runenv = {
    "docker": "ebelter/linux-tk:latest",
    "cpu": 1,
    "memory": 4,
    "disks": 20,
  }

  RunEnv bwa_runenv = {
    "docker": bwa_docker,
    "cpu": 6,
    "memory": 36,
    "disks": 20,
  }

  RunEnv samtools_runenv = {
    "docker": samtools_docker,
    "cpu": 1,
    "memory": 4,
    "disks": 20,
  }

  RunEnv freebayes_renenv = {
    "docker": freebayes_docker,
    "cpu": 1,
    "memory": 4,
    "disks": 20,
  }

  RunEnv gatk_renenv = {
    "docker": gatk_docker,
    "cpu": 4,
    "memory": 24,
    "disks": 20,
  }

  RunEnv bedtools_runenv = {
    "docker": bedtools_docker,
    "cpu": 1,
    "memory": 4,
    "disks": 20,
  }

  RunEnv abra2_renenv = {
    "docker": abra2_docker,
    "cpu": 2,
    "memory": 20,
    "disks": 20,
  }

  RunEnv dv_runenv = {
    "docker": deepvariant_docker,
    "cpu": 20,
    "memory": 96,
    "disks": 20,
  }

  call idx.run_untar_idx as reference { input:
    idx=idx,
    runenv=idx_runenv,
  }

  call align.run_bwa_mem as align { input:
    sample=sample,
    library=sample+"-lib1",
    fastqs=fastqs,
    reference=reference.path,
    runenv=bwa_runenv,
  }

  call samtools.sort as samtools_sort { input:
    bam=align.bam,
    runenv=samtools_runenv,
  } 

  if ( realign_bam ) {
    call freebayes.run_left_shift_bam as left_shift { input:
      in_bam_file=samtools_sort.sorted_bam,
      in_reference_file=reference.fasta,
      in_reference_index_file=reference.fai,
      runenv=freebayes_renenv,
    }

    call samtools.index as left_shift_index { input:
      bam=left_shift.output_bam_file,
      runenv=samtools_runenv,
    } 

    call realigner_target_creator.run_realigner_target_creator as target_creator { input:
      in_bam_file=left_shift.output_bam_file,
      in_bam_index_file=left_shift_index.bai,
      in_reference_file=reference.fasta,
      in_reference_index_file=reference.fai,
      in_reference_dict_file=reference.dict,
      runenv=gatk_renenv,
    }

    call bedtools.run_slop as expand_targets { input:
      bed_file=target_creator.output_target_bed_file,
      reference_fai=reference.fai,
      bases=targets_expansion_bases,
      runenv=bedtools_runenv,
    }

    call abra2.run_realigner as realign { input:
      in_bam_file=left_shift.output_bam_file,
      in_bam_index_file=left_shift_index.bai,
      in_target_bed_file=expand_targets.slopped_bed_file,
      in_reference_file=reference.fasta,
      in_reference_index_file=reference.fai,
      runenv=abra2_renenv,
    }
  }

  if ( !realign_bam ) {
    call samtools.index as align_index { input:
      bam=align.bam,
      runenv=samtools_runenv,
    } 
  }

  File final_bam = select_first([realign.indel_realigned_bam, samtools_sort.sorted_bam])
  File final_bai = select_first([realign.indel_realigned_bam_index, align_index.bai])

 call deepvariant.run_deepvariant as dv { input:
  sample=sample,
    bam=final_bam,
    bai=final_bai,
    reference_path=reference.path,
    runenv=dv_runenv,
  }

  call samtools.stat as samtools_stat { input:
    bam=final_bam,
    runenv=samtools_runenv,
  } 

  output {
    File bam = final_bam
    File bai = final_bai
    File bam_stats = samtools_stat.stats
    File vcf = dv.vcf
    File vcf_tvi = dv.vcf_tbi
    File dv_report = dv.report
  }
}
