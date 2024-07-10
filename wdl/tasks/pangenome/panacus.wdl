version development

import "../../structs/runenv.wdl"

task run_panacus_hist {
  input {
     File gfa_gz
     RunEnv runenv
  }

  String gfa_bn = basename(gfa_gz)
  command <<<
    ln ~{gfa_gz} ~{gfa_bn}
    run_panacus_hist ~{gfa_bn}
  >>>

  output {
    File hist = glob("*.hist")[0]
    File png = glob("*.hist.png")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "${runenv.memory} GB"
    #disks: runenv.disks
  }
}
