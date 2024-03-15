version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/abra2.wdl"
import "wdl/tasks/bwa/idx.wdl"
import "wdl/tasks/deepvariant.wdl"
import "wdl/tasks/freebayes.wdl"
import "wdl/tasks/gatk/realigner_target_creator.wdl"
import "wdl/tasks/samtools.wdl"

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
    Int rtc_expansion_bases = 160
    # dockers
    String abra2_docker = "mgibio/abra2:v2.24-focal"
    String deepvariant_docker = "google/deepvariant:1.6.0"
    String freebayes_docker = "mgibio/freebayes:1.3.6-focal"
    String gatk_docker = "broadinstitute/gatk:4.3.0.0"
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

  call realigner_target_creator.run_realigner_target_creator as targets { input:
    in_bam_file=left_shift.output_bam_file,
    in_bam_index_file=left_shift_index.bai,
    in_reference_file=reference.fasta,
    in_reference_index_file=reference.fai,
    in_reference_dict_file=reference.dict,
    in_expansion_bases=rtc_expansion_bases,
    runenv=gatk_renenv,
  } 

  call abra2.run_realigner as realign { input:
    in_bam_file=left_shift.output_bam_file,
    in_bam_index_file=left_shift_index.bai,
    in_target_bed_file=targets.output_target_bed_file,
    in_reference_file=reference.fasta,
    in_reference_index_file=reference.fai,
    runenv=abra2_renenv,
  } 

  call deepvariant.run_deepvariant as dv { input:
    sample=sample,
    bam=realign.indel_realigned_bam,
    bai=realign.indel_realigned_bam_index,
    reference=reference.path,
    runenv=dv_runenv,
  }

  output {
    File vcf = dv.vcf
    File vcf_tbi = dv.vcf_tbi
    File vcf_report = dv.report
  }
}
