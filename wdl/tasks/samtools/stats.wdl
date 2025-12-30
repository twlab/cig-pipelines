version development

import "../../structs/runenv.wdl"

task run_stat {
  input {
    File bam
    RunEnv runenv
  }

  String stat_fn = "~{basename(bam)}.stat"
  command <<<
    samtools stat -@ ~{runenv.cpu} ~{bam} > ~{stat_fn}
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
