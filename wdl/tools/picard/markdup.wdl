version development

# Picard Mark Dups

import "wdl/structs/runenv.wdl"
import "wdl/tasks/picard/markdup.wdl"

workflow picard_markdup {
    input {
        String sample
        File bam
        String docker = "ebelter/picard:2.27.4"
        Int cpu = 4
        Int memory = 20
        Int disks = 20
    }

    RunEnv runenv = {
      "docker": docker,
      "cpu": cpu,
      "memory": memory,
      "disks": disks,
    }

    call markdup.run_markdup { input:
        sample=sample,
        bam=bam,
        runenv=runenv,
    }

    output {
        File dedup_bam = run_markdup.dedup_bam
        File metrics = run_markdup.metrics
    }
}
