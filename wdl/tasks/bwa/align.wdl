version development

import "../../structs/runenv.wdl"

task run_bwa_mem {
  input {
    String sample
    String library
    Array[File] fastqs
    Array[File] idx_files
    RunEnv runenv
  }

  String bam = "~{sample}.bam"
  Int bwa_cpu = runenv.cpu - 1
  command <<<
    mkdir ref
    for f in ~{sep=' ' idx_files}; do
      ln ${f} ref/
    done
    reference_fasta=$(find ref/ -maxdepth 1 -type f -name \*.fasta)
    bwa mem \
      -t ~{bwa_cpu} \
      -K 320000000 \
      -R '@RG\tID:~{library}\tLB:~{library}\tSM:~{sample}\tPL:illumina' \
      $reference_fasta \
      ~{fastqs[0]} \
      ~{default="" fastqs[1]} | \
      samtools view -hbS - > ~{bam}
  >>>

  runtime {
    docker: runenv.docker
    cpu : runenv.cpu
    memory : "~{runenv.memory} GB"
  }

  output {
    File bam = "~{bam}"
  }
}

task run_bwamem_with_sort {
  input {
    String sample
    String library
    Array[File] fastqs
    Array[File] idx_files
    RunEnv runenv
  }

  String bam = "~{sample}.bam"
  Int bwa_cpu = runenv.cpu - 2
  Int sort_cpu = runenv.cpu - bwa_cpu
  command <<<
    mkdir ref
    for f in ~{sep=' ' idx_files}; do
      ln ${f} ref/
    done
    reference_fasta=$(find ref/ -maxdepth 1 -type f -name \*.fasta)
    bwa mem \
      -t ~{bwa_cpu} \
      -K 320000000 \
      -R '@RG\tID:~{library}\tLB:~{library}\tSM:~{sample}\tPL:illumina' \
      $reference_fasta \
      ~{fastqs[0]} \
      ~{default="" fastqs[1]} | \
      samtools sort -@ ~{sort_cpu} -O bam -o ~{bam} -
  >>>

  runtime {
    docker: runenv.docker
    cpu : runenv.cpu
    memory : "~{runenv.memory} GB"
  }

  output {
    File bam = "~{bam}"
  }
}
