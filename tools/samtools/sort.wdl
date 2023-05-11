version development

import "../../structs/runenv.wdl"
import "../../tasks/samtools.wdl"

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

    call samtools.sort as sort { input:
        bam=bam,
        runenv=runenv,
    }

    output {
        File sorted_bam = sort.sorted_bam
    }
}
