version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/vcallers/deepvariant.wdl"

workflow deepvariant {
  input {
    File sample
    File bam
    File bai
    File reference_fasta
    File reference_fai
    File reference_dict
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

  call deepvariant.run_deepvariant as dv { input:
    sample=sample,
    bam=bam,
    bai=bai,
    ref_fasta=reference_fasta,
    ref_fai=reference_fai,
    ref_dict=reference_dict,
    runenv=runenv,
  }
}
