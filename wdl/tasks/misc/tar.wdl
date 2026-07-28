version development

import "../../structs/runenv.wdl"

task run_tar {
  input {
    String output_file
    Array[File] files
    RunEnv runenv
  }

  command <<<
    set -e
    mkdir working
    cd working
    for f in ~{sep=" " files}; do
      ln ${f} .
    done
    tar cvvf ../~{output_fn} *
  >>>

  output {
    File tar_file = output_file
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
    disks: "local-disk ~{runenv.disks} SSD"
  }
}
