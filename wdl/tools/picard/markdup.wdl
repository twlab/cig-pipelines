version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/picard/markdup.wdl"

workflow picard_markdup {
  input {
    File bam
    String params
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

  call markdup.run_markdup { input:
    bam=bam,
    params=params,
    runenv=runenv,
  }

  output {
    File dedup_bam = run_markdup.dedup_bam
    File metrics = run_markdup.metrics
  }
}
