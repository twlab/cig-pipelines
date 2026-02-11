version development

import "../../structs/runenv.wdl"

task run_extract_reads {
  input {
    File sam
    File? ref
    String output_fof
    RunEnv runenv
  }

  command <<<
    samtools view -@ ~{runenv.cpu - 1} ~{if (defined(ref)) then "--reference ~{ref}" else ""} ~{sam} | awk '{print $1}' > ~{output_fof}
  >>>

  output {
    File fof = output_fof
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
  }
}
