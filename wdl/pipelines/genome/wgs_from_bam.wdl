version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/abra2.wdl"
import "wdl/tasks/bed/bedtools.wdl"
import "wdl/tasks/bwa/idx.wdl"
import "wdl/tasks/freebayes.wdl"
import "wdl/tasks/gatk/realigner_target_creator.wdl"
import "wdl/tasks/samtools.wdl"
import "wdl/tasks/vcallers/deepvariant.wdl"

workflow genome_wgs_from_bam {
  meta {
      author: "Eddie Belter"
      version: "1.2"
      description: "Starting with a BAM, Left/Realign Indels then Call Variants with DeepVariant"
  }

  input {
    String sample
    File sorted_bam                    # must be coordiante sorted
    File idx                           # tarred BWA index
    Int targets_expansion_bases = 160
    # dockers
    String abra2_docker = "mgibio/abra2:v2.24-focal"
    String bwa_docker = "mgibio/bwa:0.7.17"
    String bedtools_docker = "biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1"
    String deepvariant_docker = "google/deepvariant:1.5.0"
    String freebayes_docker = "mgibio/freebayes:1.3.6-focal"
    String gatk_docker = "broadinstitute/gatk3@sha256:5ecb139965b86daa9aa85bc531937415d9e98fa8a6b331cb2b05168ac29bc76b" #"broadinstitute/gatk:4.3.0.0"
    String samtools_docker = "mgibio/samtools:1.15.1-buster"
  }

  # RunEnvs in order of usage
  RunEnv idx_runenv = {
    "docker": "mgibio/linux-tk:latest",
    "cpu": 1,
    "memory": 4,
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

  call freebayes.run_left_shift_bam as left_shift { input:
    in_bam_file=sorted_bam,
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

  call deepvariant.run_deepvariant as dv { input:
    sample=sample,
    bam=realign.indel_realigned_bam,
    bai=realign.indel_realigned_bam_index,
    ref_fasta=reference.fasta,
    ref_fai=reference.fai,
    ref_dict=reference.dict,
    runenv=dv_runenv,
  }

  call samtools.stat as samtools_stat { input:
    bam=realign.indel_realigned_bam,
    runenv=samtools_runenv,
  } 

  output {
    File bam = realign.indel_realigned_bam
    File bai = realign.indel_realigned_bam_index
    File bam_stats = samtools_stat.stats
    File vcf = dv.vcf
    File vcf_tvi = dv.vcf_tbi
    File dv_report = dv.report
  }
}
