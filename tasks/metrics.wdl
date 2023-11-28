version development

import "../structs/runenv.wdl"

task run_seqlendist { # requires CIG metrics command
  input {
    Array[File] seqfiles
    Array[String] labels
    Array[String] reports
    String bins = "lr"
    String out
    RunEnv runenv
  }

  command <<<
    metrics seqlendist -o ~{out} -l ~{sep="," labels} -r ~{sep=" -r " reports} -b ~{bins} ~{sep=" " seqfiles}
  >>>

  output {
    Array[File] reports = glob("~{out}.*")
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    #disks : select_first([runenv.disks,"local-disk 100 SSD"])
  }
}
