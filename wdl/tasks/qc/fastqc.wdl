version development

import "../structs/runenv.wdl"

task run_fastqc {
  input {
    Array[File] seqfiles
    String format = "fastq" # sam bam
    RunEnv runenv
  }

  String output_dir = "fastqc"
  command <<<
    mkdir ~{output_dir}
    fastqc ~{sep=" " seqfiles} -f ~{format} -o ~{output_dir}
  >>>

  output {
    Array[File] output_files = glob("~{output_dir}/*.zip")
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks: runenv.disks
  }
}
