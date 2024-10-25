version development

import "../../structs/runenv.wdl"

task generate_dist {
  input {
     File gbz
     RunEnv runenv
  }

  String bn = basename(gbz, ".gbz")
  String dist = bn + ".dist"
  Int threads = runenv.cpu
  command <<<
    vg index --progress --threads ~{threads} --temp-dir ./tmp/ -j ~{dist} ~{gbz}
  >>>

  output {
    File dist = glob("~{dist}")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
    #disks: runenv.disks
  }
}
