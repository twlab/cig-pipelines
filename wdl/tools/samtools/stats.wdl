version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/samtools/stats.wdl"

workflow samtools_stats {
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

  call stats.run_stats { input:
    sam_file=sam_file,
    reference=reference,
    runenv=runenv
  }

  output {
    File stats = run_stats.stats
  }
}
