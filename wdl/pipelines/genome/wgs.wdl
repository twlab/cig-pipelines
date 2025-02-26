version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/abra2.wdl"
import "wdl/tasks/bed/bedtools.wdl"
import "wdl/tasks/bwa/align.wdl"
import "wdl/tasks/bwa/idx.wdl"
import "wdl/tasks/picard/markdup.wdl"
import "wdl/tasks/qc/fastqc.wdl"
import "wdl/tasks/freebayes.wdl"
import "wdl/tasks/gatk/realigner_target_creator.wdl"
import "wdl/tasks/samtools.wdl"
import "wdl/tasks/vcallers/deepvariant.wdl"

workflow genome_wgs {
  meta {
      author: "Eddie Belter"
      version: "1.2"
      description: "Align with BWA, Left/Realign Indels then Call Variants with DeepVariant"
  }

  input {
    String sample
    Array[File] fastqs
    File idx # tarred BWA index
    Boolean generate_fastqc = false
    Boolean realign_bam = true
    Int targets_expansion_bases = 160
    # dockers and resources
    String abra2_docker
    Int abra2_cpu
    Int abra2_memory
    String bwa_docker
    Int bwa_cpu
    Int bwa_memory
    String bedtools_docker
    Int bedtools_cpu
    Int bedtools_memory
    String deepvariant_docker
    Int deepvariant_cpu
    Int deepvariant_memory
    String? fastqc_docker
    Int? fastqc_cpu
    Int? fastqc_memory
    String freebayes_docker
    Int freebayes_cpu
    Int freebayes_memory
    String gatk_docker
    Int gatk_cpu
    Int gatk_memory
    String markdup_docker
    Int markdup_cpu
    Int markdup_memory
    String samtools_docker
    Int samtools_cpu
    Int samtools_memory
    String utils_docker
    Int utils_cpu
    Int utils_memory
  }

  # RunEnvs in order of usage
  RunEnv utils_runenv = {
    "docker": utils_docker,
    "cpu": utils_cpu,
    "memory": utils_memory,
    "disks": 20,
  }

  RunEnv bwa_runenv = {
    "docker": bwa_docker,
    "cpu": bwa_cpu,
    "memory": bwa_memory,
    "disks": 20,
  }

  RunEnv samtools_runenv = {
    "docker": samtools_docker,
    "cpu": samtools_cpu,
    "memory": samtools_memory,
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

  RunEnv markduper_runenv = {
    "docker": markdup_docker,
    "cpu": markdup_cpu,
    "memory": markdup_memory,
    "disks": 20,
  }

  RunEnv dv_runenv = {
    "docker": deepvariant_docker,
    "cpu": deepvariant_cpu,
    "memory": deepvariant_memory,
    "disks": 20,
  }

  if ( generate_fastqc ) {
    RunEnv runenv_fastqc = {
      "docker": fastqc_docker,
      "cpu": fastqc_cpu,
      "memory": fastqc_memory,
      "disks": 20,
    }
    call fastqc.run_fastqc { input:
      seqfiles=fastqs,
      runenv=runenv_fastqc,
    }
  }

  call idx.run_untar_idx as reference { input:
    idx=idx,
    runenv=utils_runenv,
  }

  call align.run_bwa_mem as align { input:
    sample=sample,
    library=sample+"-lib1",
    rg_id=sample+"-lib1",
    platform_unit=sample+"-lib1",
    fastqs=fastqs,
    idx_files=[reference.fasta, reference.amb, reference.ann, reference.bwt, reference.pac, reference.sa],
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

  call markdup.run_markdup as picard_markdup { input:
    bam=select_first([realign.indel_realigned_bam, samtools_sort.sorted_bam]), # select realign first if true
    runenv=markduper_runenv,
  }

  call samtools.index as samtools_index { input:
    bam=picard_markdup.dedup_bam,
    runenv=samtools_runenv,
  }

  call deepvariant.run_deepvariant as dv { input:
    sample=sample,
    bam=picard_markdup.dedup_bam,
    bai=samtools_index.bai,
    ref_fasta=reference.fasta,
    ref_fai=reference.fai,
    ref_dict=reference.dict,
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
