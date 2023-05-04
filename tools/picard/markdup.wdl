version development

# Picard Mark Dups

import "wdl/structs/runenv.wdl"
import "wdl/tasks/picard/markdup.wdl"

workflow markdup {
    input {
        String name
        File bam
        String docker
    }

    RunEnv runenv = {
      "docker": docker,
      "cpu": 4,
      "memory": 20,
    }

    call picard.run_markdup { input:
        name=name,
        bam=bam,
        runenv=runenv,
    }

    output {
        File bam = run_markdup.dedup_bam
        File metrics = run_markdup.metrics
    }
}
