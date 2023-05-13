version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/samtools.wdl"

workflow index {
    input {
        File bam
        String docker = "ebelter/samtools:1.15.1"
        Int cpu = 1
        Int memory = 4
        Int disks = 20
    }

    RunEnv runenv = {
      "docker": docker,
      "cpu": cpu,
      "memory": memory,
      "disks": disks,
    }

    call samtools.index as index { input:
        bam=bam,
        runenv=runenv
    }

    output {
        File bai = index.bai
    }
}
