version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/abra2.wdl"
import "wdl/tasks/bed/bedtools.wdl"
import "wdl/tasks/bwa/idx.wdl"
import "wdl/tasks/freebayes.wdl"
import "wdl/tasks/gatk/realigner_target_creator.wdl"
import "wdl/tasks/samtools.wdl"
import "wdl/tasks/vcallers/deepvariant.wdl"

workflow realign_and_dv {
  meta {
    author: "Eddie Belter"
    version: "0.1"
    description: "Realign Reads and Run DeepVariant Pipleine"
  }

  input {
    String sample
    File bam
    File bai
    File idx
    Int targets_expansion_bases = 160
    # dockers
    String abra2_docker = "mgibio/abra2:v2.24-focal"
    String bedtools_docker = "biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1"
    String deepvariant_docker = "google/deepvariant:1.6.0"
    String freebayes_docker = "mgibio/freebayes:1.3.6-focal"
    String gatk_docker = "broadinstitute/gatk3@sha256:5ecb139965b86daa9aa85bc531937415d9e98fa8a6b331cb2b05168ac29bc76b" #"broadinstitute/gatk:4.3.0.0"
    String samtools_docker = "mgibio/samtools:1.15.1-buster"
  }

  # RunEnvs in order of usage
  RunEnv runenv_idx = {
    "docker": "ebelter/linux-tk:latest",
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

  RunEnv samtools_runenv = {
    "docker": samtools_docker, 
    "cpu": 1,
    "memory": 4,
    "disks": 20,
  }

  RunEnv gatk_renenv = {
    "docker": gatk_docker,
    "cpu": 1,
    "memory": 4,
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
    "cpu": 1,
    "memory": 4,
    "disks": 20,
  }

  RunEnv dv_runenv = {
    "docker": deepvariant_docker,
    "cpu": 9,
    "memory": 48,
    "disks": 20,
  }

  # Calls
  call idx.run_untar_idx as reference { input:
    idx=idx,
    runenv=runenv_idx,
  }

  call freebayes.run_left_shift_bam as left_shift { input:
    in_bam_file=bam,
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
    reference_path=reference.path,
    runenv=dv_runenv,
  }

  output {
    File vcf = dv.vcf
    File vcf_tbi = dv.vcf_tbi
    File vcf_report = dv.report
  }
}
