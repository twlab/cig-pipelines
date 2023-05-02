version development

# Picard Metrics for RNA Seq with helper tools
# these tasks require picard

import "../../structs/runenv.wdl"

task create_ribosomal_intervals {
    input {
        File alignments
        File annotation
        RunEnv runenv
    }

    String intervals = "rRNA.intervals_list"
    command <<<
        create_ribosomal_intervals \
            ~{alignments} \
            ~{annotation} \
            ~{intervals}
    >>>

    output {
        File intervals = "${intervals}"
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: "~{select_first([runenv.memory,'8'])} GB"
    }
}

task create_refflat {
    input {
        File annotation
        RunEnv runenv
    }

    command <<<
        create_refflat \
            ~{annotation} \
            refflat
    >>>

    output {
        File refflat = "refflat"
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: "~{select_first([runenv.memory,'8'])} GB"
    }
}

task run_collect_rnaseq_metrics {
    input {
        File alignments
        File ribosomal_intervals
        File refflat
        RunEnv runenv
    }

    String metrics = basename(alignments) + ".metrics"
    command <<<
        java -jar /usr/picard/picard.jar CollectRnaSeqMetrics \
            --INPUT ~{alignments} \
            --STRAND NONE \
            --REF_FLAT ~{refflat} \
            --RIBOSOMAL_INTERVALS ~{ribosomal_intervals} \
            --OUTPUT ~{metrics}
    >>>

    output {
        File metrics = "${metrics}"
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: "~{select_first([runenv.memory,'8'])} GB"
    }
}
