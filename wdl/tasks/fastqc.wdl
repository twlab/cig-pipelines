version development

import "../structs/runenv.wdl"

task run_fastqc {
  input {
    Array[File] seqfiles
    String format = "fastq" # sam bam
    String output_dir = "fastqc"
    RunEnv runenv
  }

  command <<<
    mkdir ~{output_dir}
    fastqc ~{sep=" " seqfiles} -f ~{format} -o ~{output_dir}
  >>>

  output {
    Directory output_dir = "~{output_dir}"
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks: runenv.disks
  }
}
