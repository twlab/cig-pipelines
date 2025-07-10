version development

import "../../structs/runenv.wdl"

task run_tar {
  input {
    String name
    Array[File] files
    RunEnv runenv
  }

  command <<<
    mkdir working
    cd working
    for f in ~{sep=" " files}; do
      ln ${f} .
    done
    tar cvvf ../~{name}.tar *
  >>>

  output {
    File tar_file = glob("~{name}.tar")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
    disks: "local-disk ~{runenv.disks} SSD"
  }
}
