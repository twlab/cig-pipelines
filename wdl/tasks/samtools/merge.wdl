version development

import "../../structs/runenv.wdl"

task run_merge_crams {
  meta {
    author: "Eddie Belter"
    version: "1.0"
    description: "Merge CRAMs."
  }

  input {
    Array[File] crams 
    File reference
    String output_cram_fn
    RunEnv runenv
  }

  command <<<
    set -e
    samtools merge -@ ~{runenv.cpu} --reference ~{reference} --write-index -O CRAM -o ~{output_cram_fn} ~{sep=' ' crams} 
  >>>

  output {
    File merged_cram = "~{output_cram_fn}"
    File merged_crai = "~{output_cram_fn}.crai"
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
  }
}


