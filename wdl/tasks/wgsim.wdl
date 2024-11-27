version development

import "../structs/runenv.wdl"

task run_wgsim {
  input {
    File fasta
    Float base_error
    Int out_distance
    Int stdev
    Int number_pairs
    Int read1_length
    Int read2_length
    Float mutation_rate
    Float fraction_indels
    Float prob_indel_extentsion
    Int seed
    RunEnv runenv
  }

  String bn = basename(fasta)
  String fq1 = "~{bn}.R1.fastq"
  String fq2 = "~{bn}.R2.fastq"
  command <<<
    wgsim -e ~{bn} -d ~{out_distance} -s ~{stdev} -N ~{number_pairs} -1 ~{read1_length} -2 ~{read2_length} -r ~{mutation_rate} -R ~{fraction_indels} -X ~{prob_indel_extentsion} -S ~{seed} ~{fasta} ~{fq1} ~{fq2}
  >>>

  output {
    Array[File] simulated_fastqs = glob("~{bn}.R*.fastq")
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    #disks : select_first([runenv.disks,"local-disk 100 SSD"])
  }
}
