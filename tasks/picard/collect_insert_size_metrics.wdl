version development

# Picard Collect Insert Size Metrics

import "../../structs/runenv.wdl"

task run_collect_insert_size_metrics {
    input {
        File alignments
        RunEnv runenv
    }

    String metrics_fn = basename(alignments) + ".insert-size.metrics.txt"
    String hist_fn = basename(alignments) +  ".insert-size.hist.pdf"
    Int javamem = runenv.memory - 2
    command <<<
        java -Xmx~{javamem}g -jar /usr/picard/picard.jar CollectInsertSizeMetrics \
            --INPUT ~{alignments} \
            --MINIMUM_PCT 0.05 \
            --OUTPUT ~{metrics_fn} \
            --Histogram_FILE ~{hist_fn}
    >>>

    output {
        File metrics = "${metrics_fn}"
        File histogram = "${hist_fn}"
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: "~{runenv.memory} GB"
    }
}
