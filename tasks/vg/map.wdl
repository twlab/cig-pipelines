version development

import "../../structs/runenv.wdl"

task run_map {
  input {
     Array[File] fastqs
     File xg
     File gcsa
     String name
     RunEnv runenv
  }

  String gam = "${name}.gam"
  # From https://github.com/vgteam/vg/wiki/Basic-Operations
  # Paired end reads in paired FASTQs
  #  vg map -x index.xg -g index.gcsa -f reads1.fq -f reads2.fq > mapped.gam
  # 
  # Params
  # graph/index:
  #  -d, --base-name BASE          use BASE.xg and BASE.gcsa as the input index pair
  #  -x, --xg-name FILE            use this xg index or graph (defaults to <graph>.vg.xg)
  #  -g, --gcsa-name FILE          use this GCSA2 index (defaults to <graph>.gcsa)
  #  -1, --gbwt-name FILE          use this GBWT haplotype index (defaults to <graph>.gbwt)
  # algorithm:
  #  -t, --threads N               number of compute threads to use
  # input:
  #  -f, --fastq FILE              input fastq or (2-line format) fasta, possibly compressed, two are allowed, one for each mate
  command <<<
    vg run -t ~{runenv.cpu - 1} -x ~{xg} -g ~{gcsa} -f ~{fastqs[0]} -f ~{fastqs[1]} > ~{gam}
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
