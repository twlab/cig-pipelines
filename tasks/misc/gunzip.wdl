version development

import "../../structs/runenv.wdl"

task run_gunzip {
  input {
    File compressed
    RunEnv runenv
  }

  String uncompressed = sub(basename(compressed), ".gz$", "")
  command <<<
    gunzip -c ~{compressed} > ~{uncompressed}
  >>>

  output {
    File uncompressed = "~{uncompressed}"
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
    disks: "local-disk ~{runenv.disks} SSD"
  }
}
