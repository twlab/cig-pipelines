version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/bwa/idx.wdl"
import "wdl/tasks/vcallers/deepvariant_pg.wdl"

workflow deepvariant_pg {
  input {
    File sample
    File bam
    File bai
    File idx
    File gbz
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

  call deepvariant_pg.run_deepvariant_pg as dv { input:
    sample=sample,
    bam=bam,
    bai=bai,
    ref_fasta=reference.fasta,
    ref_fai=reference.fai,
    ref_dict=reference.dict,
    gbz=gbz,
    runenv=dv_runenv,
  }
}
