version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/samtools.wdl"

workflow samtools_sort {
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

    call samtools.sort { input:
        bam=bam,
        runenv=runenv,
    }

    output {
        File sorted_bam = sort.sorted_bam
    }
}
