version development

import "../../structs/runenv.wdl"

task generate_haplotypes_for_hapsamp {
  input {
     File gbz
     File dist
     File r_index
     RunEnv runenv
  }

  String bn = basename(gbz, ".gbz")
  String hap = bn + ".hap"
  Int threads = runenv.cpu
  command <<<
    ln -s ~{gbz} ./
    ln -s ~{dist} ./
    ln -s ~{r_index} ./
    vg haplotypes --verbosity 2 --threads ~{threads} -r ~{r_index} -H ~{hap} ~{basename(gbz)}
  >>>

  output {
    File hap = glob("~{hap}")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
    #disks: runenv.disks
  }
}

task apply_haplotypes_to_graph {
  input {
    Array[File] fastqs
    File gbz
    File hap
    File kmers
    RunEnv runenv
  }

  Int threads = runenv.cpu
  String output_gbz = basename(gbz, ".gbz") + ".sampled.gbz"
  command <<<
    vg haplotypes --verbosity 2 --threads ~{threads} --include-reference --diploid-sampling -i ~{hap} -k ~{kmers} -g ~{output_gbz} ~{gbz}
  >>>

  output {
    File output_gbz = glob("~{output_gbz}")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "${runenv.memory} GB"
    #disks: runenv.disks
  }
}

task generater_kmers_with_kmc {
  input {
    String sample
    Array[File] fastqs
    RunEnv runenv
  }

  command <<<
    mkdir ./tmp/
    kmc -t~{runenv.cpu} -k29 -m128 -okff @~{write_lines(fastqs)} ~{sample} ./tmp/
  >>>

  output {
    File kmers = glob("*.kff")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "${runenv.memory} GB"
    #disks: runenv.disks
  }
}
