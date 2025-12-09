version development

import "../../structs/runenv.wdl"

task run_phase_crams {
  meta {
    author: "Eddie Belter"
    version: "1.0"
    description: "Merge CRAMs."
  }

  input {
    Array[File] crams 
    File reference
    File phased_reads_fof
    RunEnv runenv
  }

  command <<<
    set -e
    for cram in ~{sep=' ' crams}; do
      ln "${cram}" ./
      cram_bn=$(basename "${cram}")
      echo "${cram_bn}" >> crams.fof
    done
    /usr/local/bin/phase-reads -c crams.fof -r ~{reference} -p ~{phased_reads_fof}
  >>>

  output {
    Array[File] phased_crams = glob("*.phased.cram")
    Array[File] yamls = glob("*.yaml")
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
  }
}


