version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/pangenome/panacus.wdl"

workflow panacus {
  input {
    File gfa_gz
    String docker = "mgibio/panacus:0.2.3-buster"
    Int cpu = 1
    Int memory = 24
  }

  RunEnv runenv_panacus = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
    "disks": 20,
  }

  call panacus.run_panacus_hist { input:
    gfa_gz=gfa_gz,
    runenv=runenv_panacus,
  }
}
