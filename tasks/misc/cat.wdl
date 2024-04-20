version development

import "../../structs/runenv.wdl"

task cat {
  input {
    Array[File] files
    String out = "concat"
    RunEnv runenv
  }

  command <<<
    cat ~{sep=" " files} > ~{out}
  >>>

  output {
    File concatenated_file = out
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
    #disks: runenv.disks
  }
}

task run_zcat {
  input {
    Array[File] files
    String out = "concat"
    RunEnv runenv
  }

  command <<<
    zcat ~{sep=" " files} > ~{out}
  >>>

  output {
    File concatenated_file = out
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
    #disks: runenv.disks
  }
}
