version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/bwa/idx.wdl"
import "wdl/tasks/vcallers/deepvariant.wdl"

workflow deepvariant {
  input {
    File sample
    File bam
    File bai
    File idx
    String deepvariant_docker
    Int deepvariant_cpu
    Int deepvariant_memory
    String utils_docker
    Int utils_cpu
    Int utils_memory
  }

  RunEnv utils_runenv = {
    "docker": utils_docker,
    "cpu": utils_cpu,
    "memory": utils_memory,
    "disks": 20,
  }

  RunEnv dv_runenv = {
    "docker": deepvariant_docker,
    "cpu": deepvariant_cpu,
    "memory": deepvariant_memory,
    "disks": 20,
  }

  call idx.run_untar_idx as reference { input:
    idx=idx,
    runenv=utils_runenv,
  }

  call deepvariant.run_deepvariant as dv { input:
    sample=sample,
    bam=bam,
    bai=bai,
    ref_fasta=reference.fasta,
    ref_fai=reference.fai,
    ref_dict=reference.dict,
    runenv=dv_runenv,
  }
}
