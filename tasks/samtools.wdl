version development

# these tasks require samtools

import "../structs/runenv.wdl"

task faidx {
    input {
        File fasta
        RunEnv runenv
    }

    String bn = basename(fasta)
    command {
        samtools faidx ${fasta} --fai-idx ${bn + ".fai"} --gzi-idx ${bn + ".gzi"}
    }

    output {
        File fai = "${bn}.fai"
				File gzi = "${bn}.gzi"
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: runenv.memory + " GB"
        #disks : select_first([runenv.disks,"local-disk 100 SSD"])
    }
}

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

task merge_bams {
    input {
        String sample
        Array[File] bams
        RunEnv runenv
    }

    String output_bam = "~{sample}.bam"
    Int merge_cpu = if runenv.cpu > 1 then runenv.cpu - 1 else 1
    command <<<
        samtools merge -n -r -@ ~{merge_cpu} -o ~{output_bam} ~{sep=" " bams}
    >>>

    output {
        File merged_bam = output_bam
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
        ln ~{bam} ~{basename(bam)}
        samtools index -b -@~{samtools_cpu} ~{basename(bam)}
    >>>

    output {
        File bai = "~{basename(bam)}.bai"
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
    samtools sort -O ~{output_fmt} -o ~{output_bam} ~{bam}
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
        File stats = "${stat_fn}"
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: "~{runenv.memory} GB"
    }
}
