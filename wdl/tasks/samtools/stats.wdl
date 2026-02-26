version development

import "../../structs/runenv.wdl"

task run_stats {
  input {
    File sam_file
    File? reference
    RunEnv runenv
  }

  String stat_fn = "~{basename(sam_file)}.stat"
  command <<<
    samtools stat -@ ~{runenv.cpu}  ~{if (defined(reference)) then "--reference ~{reference}" else ""} ~{sam_file} > ~{stat_fn}
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
