version development

# these tasks require samtools

import "../structs/runenv.wdl"

task merge_sams {
    input {
        String prefix
        File header
        Array[File] sams
        RunEnv runenv
    }

    String output_sam = prefix+".sam"
    command <<<
       cp ~{header} ~{output_sam}
       for sam in ~{sep(' ', sams)}
       do
            echo Appending ${sam}
            cat ${sam} >> ~{output_sam}
       done
    >>>

    output {
        File sam = output_sam
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: "~{runenv.memory} GB"
    }
}

task index {
    input {
        File bam
        RunEnv runenv
    }

    Int samtools_cpu = runenv.cpu
    command <<<
        set -x
        samtools index -b -@~{samtools_cpu} ~{bam}
    >>>

    output {
        File bai = "~{bam}.bai"
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: "~{runenv.memory} GB"
    }
}

task sort {
    input {
        File bam
        String output_fmt = "BAM"
        RunEnv runenv
    }

    String output_bam = basename(bam)
    command <<<
        samtools sort ~{bam} -O ~{output_fmt} -o ~{output_bam}
    >>>

    output {
        File sorted_bam = output_bam
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: "~{runenv.memory} GB"
    }
}

task stat {
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
        memory: "~{runenv.memory} GB"
    }
}
