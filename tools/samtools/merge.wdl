version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/samtools.wdl"

workflow samtools_merge {
    input {
        String name
        Array[File] bams
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

    call samtools.merge_bams as merge { input:
        name=name,
        bams=bams,
        runenv=runenv,
    }

    output {
        File merged_bam = merge.merged_bam
    }
}
