version development

import "../../structs/runenv.wdl"

task run_fastp {
  input {
    Array[File] fastqs # must be paired
    String? params
    RunEnv runenv
  }

  String trimmed_fq1 = "trimmed.R1.fastq.gz"
  String trimmed_fq2 = "trimmed.R2.fastq.gz"
  command <<<
    fastp \
      ~{params} \
      -i ~{fastqs[0]} -I ~{fastqs[1]} \
      -o ~{trimmed_fq1} -O ~{trimmed_fq2}
  >>>

  output {
    Array[File] trimmed_fastqs = glob("trimmed.*.fastq.gz")
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    #disks : select_first([runenv.disks,"local-disk 100 SSD"])
  }
}
