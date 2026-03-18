version development

import "../../structs/runenv.wdl"

task run_align {
  input {
    File query
    File target
    String output_fn = "paf"
    String params = ""
    RunEnv runenv
  }

  command <<<
    minimap2 -t ~{runenv.cpu} -o ~{output_fn} ~{params} ~{target} ~{query}
  >>>

  output {
    File alignments = output_fn
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    #disks : select_first([runenv.disks,"local-disk 100 SSD"])
  }
}
