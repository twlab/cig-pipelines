version development

import "../../structs/runenv.wdl"

task run_index {
  input {
    File bam
    RunEnv runenv
  }

  String bai = "~{basename(bam)}.bai"
  command <<<
    ln ~{bam} ~{basename(bam)}
    samtools index -b -@~{runenv.cpu} -o ~{bai} ~{basename(bam)}
  >>>

  output {
    File bai = glob(bai)[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
  }
}
