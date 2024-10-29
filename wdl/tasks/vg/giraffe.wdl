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

task run_giraffe_haplotype_mode {
  input {
    String sample
    Array[File] fastqs
    File gbz
    File haplotypes
    File kmers
    RunEnv runenv
  }

  String gam = "${sample}.gam"
  # input options:
  #  -f  fastqs
  # basic options
  #  -Z  use this GBZ file (GBWT index + GBWTGraph)
  #  -p  show progress
  #  -o  output the alignments in NAME format [gam]
  #  -t  number of mapping threads to use
  # haplotype sampling:
  #  --haplotype-name FILE         sample from haplotype information in FILE
  #  --kff-name FILE               sample according to kmer counts in FILE
  command <<<
    vg giraffe -p -t ~{runenv.cpu - 1} --kff-name ~{kmers} --haplotype-name ~{haplotypes} -Z ~{gbz} -f ~{fastqs[0]} -f ~{fastqs[1]} -o gam > ~{gam}
  >>>

  output {
    File gam = glob(gam)[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "${runenv.memory} GB"
    #disks: runenv.disks
  }
}
