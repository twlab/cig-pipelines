version development

# Picard Mark Dups
# these tasks require picard

import "../../structs/runenv.wdl"

task run_markdup {
    input {
        String name
        File bam
        RunEnv runenv
    }

    String output_bam = "${name}.dedup.bam"
    String output_metrics = "${name}.dedup.metrics"
    Int javamem = runenv.memory - 2
    command <<<
        java -Xmx~{javamem}g -jar /usr/picard/picard.jar MarkDuplicates \
            --INPUT ~{bam} \
            --OUTPUT ~{output_bam} \
            --METRICS_FILE ~{output_metrics} \
            --QUIET true \
            --VALIDATION_STRINGENCY LENIENT
    >>>

    output {
        File dedup_bam = "${output_bam}"
        File metrics = "${output_metrics}"
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: runenv.memory + " GB"
    }
}
