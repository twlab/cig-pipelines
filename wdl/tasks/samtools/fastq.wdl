version development

import "../../structs/runenv.wdl"

task run_sam_to_fastq {
  input {
    File sam
    File? reference # fasta (gz) then optional indexes - these are created by satools if not available
    String params = ""
    RunEnv runenv
  }

  String fastq = sub(basename(sam), "\\..*", ".fastq")
  command <<<
    samtools fastq -@ ~{runenv.cpu} ~{params} -o ~{fastq} ~{if (defined(reference)) then "--reference ~{reference}" else ""} ~{sam}
  >>>

  output {
    File fastq = "~{fastq}"
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
  }
}
