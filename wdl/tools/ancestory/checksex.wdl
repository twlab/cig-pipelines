version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/checksex.wdl"

workflow checksex {
  meta {
    author: "Eddie Belter"
    version: "1.0"
  }

  input {
    File input_vcf
    File input_vcf_tbi
    String docker = "mgibio/smahtqc:v1.0-buster"
    Int cpu = 1
    Int memory = 8
  }

  RunEnv runenv = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
    "disks": 20,
  }

  call checksex.run_checksex { input:
    input_vcf=input_vcf,
    input_vcf_tbi=input_vcf_tbi,
    runenv=runenv,
  }

  output {
    File output_file = run_checksex.output_file
  }
}
