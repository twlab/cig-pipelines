version development

import "../../structs/runenv.wdl"

task run_giraffe {
  input {
     Array[File] fastqs
     File gbz
     File min
     File dist
     String sample
     RunEnv runenv
  }

  String gam = "${sample}.gam"
  # -f  fastqs
  # -m  minimizer index
  # -d  distance index
  # -Z  use this GBZ file (GBWT index + GBWTGraph)
  # -p  show progress
  # -o  output the alignments in NAME format
  # -t  number of mapping threads to use
  command <<<
    vg giraffe -t ~{runenv.cpu - 1} -m ~{min} -d ~{dist} -Z ~{gbz} -f ~{fastqs[0]} -f ~{fastqs[1]} -o gam > ~{gam}
  >>>

  output {
    File gam = gam
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "${runenv.memory} GB"
    #disks: runenv.disks
  }
}
