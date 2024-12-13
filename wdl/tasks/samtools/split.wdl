version development

# these tasks require samtools

import "../../structs/runenv.wdl"

task run_split_by_chromosome {
  input {
    File fasta
    File fai
    Array[String] chrs
    RunEnv runenv
  }

  String bn = basename(fasta)
  command <<<
    for chr in ~{sep=' ' chrs}; do
      samtools faidx ~{fasta} --fai-idx ~{fai} "${chr}" -o "~{bn}.${chr}.fasta"
    done
  >>>

  output {
    Array[File] chromosome_fastas = glob("~{bn}.*.fasta")
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    #disks : select_first([runenv.disks,"local-disk 100 SSD"])
  }
}
