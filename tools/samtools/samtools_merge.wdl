version development

import "structs/runenv.wdl"

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
        cpu: runenv.cpu
        memory: "~{runenv.memory} GB"
        docker: "ebelter/samtools:1.15.1"
        disks: runenv.disks
        docker: runenv.docker
    }
}
