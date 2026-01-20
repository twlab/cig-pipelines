version development

import "../../structs/runenv.wdl"

task run_genomecov {
  meta {
    author: "Eddie Belter"
    version: "0.1"
    description: "Run Bedtools Genome Coverage"
  }

  input {
    File cram
    File ref
    RunEnv runenv
  }

  String output_file = basename(cram, ".cram") + ".bed"
  command <<<
    CRAM_REFERENCE=~{ref} bedtools genomecov -ibam ~{cram} -bga > ~{output_file}
  >>>

  output {
    File bed = output_file
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
  }
}
