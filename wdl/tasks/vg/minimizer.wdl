version development

import "../../structs/runenv.wdl"

task generate_min {
  input {
     File gbz
     File dist
     RunEnv runenv
  }

  String bn = basename(gbz, ".gbz")
  String min = bn + ".min"
  Int threads = runenv.cpu
  command <<<
    vg minimizer --progress --threads ~{threads} -o ~{min} -d ~{dist} ~{gbz} 
  >>>

  output {
    File min = glob("~{min}")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
    #disks: runenv.disks
  }
}
