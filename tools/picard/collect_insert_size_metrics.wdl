version development

# Picard Collect Insert Metrics

import "../../structs/runenv.wdl"
import "../../tasks/picard/collect_insert_size_metrics.wdl"

workflow collect_insert_size_metrics {
    input {
        File alignments
        String docker = "ebelter/picard:2.27.4"
        Int memory = 4 
        Int cpu = 20
    }

    RunEnv runenv = {
      "docker": docker,
      "cpu": cpu,
      "memory": memory,
      "disks": 20,
    }

    call collect_insert_size_metrics.run_collect_insert_size_metrics { input:
        alignments=alignments,
        runenv=runenv,
    }

    output {
        File metrics = run_collect_insert_size_metrics.metrics
        File histogram = run_collect_insert_size_metrics.histogram
    }
}
