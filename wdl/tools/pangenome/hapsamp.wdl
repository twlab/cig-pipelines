version development

import "wdl/tasks/vg/gbwt.wdl"
import "wdl/tasks/vg/haplotypes.wdl"
import "wdl/structs/runenv.wdl"

workflow pangenome_hapsamp {
  input {
    File dist
    File gbz
    File min
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

  call gbwt.generate_r_index { input:
    gbz=gbz,
    runenv=runenv,
  }

  call haplotypes.generate_haplotypes_for_hapsamp { input:
    gbz=gbz,
    dist=dist,
    r_index=generate_r_index.r_index,
    runenv=runenv,
  }
}
