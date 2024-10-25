version development

import "../../structs/runenv.wdl"

task generate_r_index {
  input {
     File gbz
     RunEnv runenv
  }

  String bn = basename(gbz, ".gbz")
  String ri = bn + ".ri"
  Int threads = runenv.cpu
  command <<<
    vg gbwt --progress --num-threads ~{threads} -r ~{ri} -Z ~{gbz}
  >>>

  output {
    File r_index = glob("~{ri}")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
    #disks: runenv.disks
  }
}
