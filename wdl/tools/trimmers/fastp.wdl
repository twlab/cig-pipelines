version development

import "wdl/tasks/trimmers/fastp.wdl"
import "wdl/structs/runenv.wdl"

workflow fastp_trimmer {
  input {
    Array[File] fastqs # must be paired
    String params = ""
    String docker
    Int cpu
    Int memory
  }

  RunEnv trimmer_runenv = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
    "disks": 20,
  }

  call fastp.run_fastp as trimmer { input:
    fastqs=fastqs,
    params=params,
    runenv=trimmer_runenv,
  }

  output {
    Array[File] trimmed_fastqs = trimmer.trimmed_fastqs
  }
}
