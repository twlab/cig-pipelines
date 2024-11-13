version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/vcallers/deepvariant.wdl"

workflow deep_variant {
  meta {
    author: "Eddie Belter"
    version: "1.2"
    description: "Call variants with Deep Variant"
  }

  input {
    String sample
    File bam
    File bai
    File ref_fasta
    File ref_fai
    File ref_dict
    String docker
    Int cpu
    Int memory
  }

  RunEnv runenv = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
    "disks": 20,
  }

  call deepvariant.NEW_run_deepvariant as dv { input:
    sample=sample,
    bam=bam,
    bai=bai,
    ref_fasta=ref_fasta,
    ref_fai=ref_fai,
    ref_dict=ref_dict,
    runenv=runenv,
  }

  output {
    File vcf = dv.vcf
    File vcf_tbi = dv.vcf_tbi
    File gvcf = dv.gvcf
    File gvcf_tbi = dv.gvcf_tbi
    File report = dv.report
  }
}
