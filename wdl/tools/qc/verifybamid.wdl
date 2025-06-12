version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/qc/contamination/verifybamid.wdl"

workflow verifybamid {
  input {
    String sample
    File bam
    File bai
    File ref_fasta
    File ref_fai
    File resource
    String docker
    Int cpu
    Int memory
  }

  RunEnv runenv_verifybamid = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
    "disks": 20,
  }

  call verifybamid.run_verifybamid { input:
    sample=sample,
    bam=bam,
    bai=bai,
    ref_fasta=ref_fasta,
    ref_fai=ref_fai,
    resource=resource,
    runenv=runenv_verifybamid,
  }
}
