version development

# Picard Mark Dups

import "../../structs/runenv.wdl"
import "../../tasks/picard/markdup.wdl"

workflow picard_markdup {
    input {
        String name
        File bam
        String docker
        Int cpu = 4
        Int memory = 20
    }

    RunEnv runenv = {
      "docker": docker,
      "cpu": cpu,
      "memory": memory,
    }

    call markdup.run_markdup { input:
        name=name,
        bam=bam,
        runenv=runenv,
    }

    output {
        File dedup_bam = run_markdup.dedup_bam
        File metrics = run_markdup.metrics
    }
}
