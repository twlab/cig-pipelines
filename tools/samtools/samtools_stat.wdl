version development

import "structs/runenv.wdl"

task run_stat {
    input {
        File bam
        RunEnv runenv
    }

    String stat_fn = "~{basename(bam)}.stat"
    command <<<
        samtools stat ~{bam} > ~{stat_fn}
    >>>

    output {
        File stat_file = "${stat_fn}"
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        disks: runenv.disks
        memory: runenv.memory + " GB"
        #disks : select_first([runenv.disks,"local-disk 100 SSD"])
    }
}
