version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/samtools/index.wdl"

workflow samtools_index {
  input {
    File sam_file
    File? reference
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

  call index.new_run_index as run_index { input:
    sam_file=sam_file,
    reference=reference,
    runenv=runenv
  }

  output {
    File idx = run_index.idx
  }
}
