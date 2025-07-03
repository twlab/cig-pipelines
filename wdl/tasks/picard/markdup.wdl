version development

import "../../structs/runenv.wdl"

task run_markdup {
    input {
        File bam
        String sort_order = "coordinate"
        RunEnv runenv
    }

    String output_bam_bn = basename(bam, ".bam")
    String output_bam = "~{output_bam_bn}.dedup.bam"
    String output_metrics = "~{output_bam_bn}.dedup.metrics"
    command <<<
        java -Xmx~{runenv.memory - 2}g -jar /usr/picard/picard.jar MarkDuplicates \
            --INPUT ~{bam} \
            --OUTPUT ~{output_bam} \
            --METRICS_FILE ~{output_metrics} \
            --QUIET true \
            --ASSUME_SORT_ORDER ~{sort_order} \
            --VALIDATION_STRINGENCY SILENT
    >>>

    output {
        File dedup_bam = "~{output_bam}"
        File metrics = "~{output_metrics}"
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: runenv.memory + " GB"
    }
}
