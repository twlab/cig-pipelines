version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/star/signals.wdl"

workflow star_signals {
  input {
    File bam
    String strandedness
    String docker
    Int cpu
    Int mem
  }

  RunEnv runenv = {
    "docker": docker,
    "cpu": cpu,
    "memory": mem,
    "disks": 20,
  }

  call signals.run_star_signals { input:
    bam=bam,
    strandedness=strandedness,
    runenv=runenv,
  }
}
