version development

import "../../structs/runenv.wdl"

task get_chromosomes {
  input {
    File fai
    RunEnv runenv
  }

  command <<<
    cut -f1 ~{fai} | tee chromosomes
  >>>

  output {
    Array[String] chromosomes = read_lines("chromosomes")
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    #disks : select_first([runenv.disks,"local-disk 100 SSD"])
  }
}
