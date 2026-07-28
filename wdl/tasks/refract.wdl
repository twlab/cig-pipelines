version development

import "../structs/runenv.wdl"

task build_idx {
  input {
    File paf
    RunEnv runenv
  }

  String output_rfx = "~{basename(paf)}.rfx"
  command <<<
    refract build-liftover-index ~{paf} ~{output_rfx}
  >>>

  output {
    File rfx = output_rfx
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
  }
}
