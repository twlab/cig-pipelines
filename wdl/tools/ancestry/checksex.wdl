version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/qc/checksex.wdl"

workflow checksex {
  meta {
    author: "Eddie Belter"
    version: "1.0"
  }

  input {
    File input_vcf
    File input_vcf_tbi
    String seq_platform
    String docker
    Int cpu
    Int memory
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
    seq_platform=seq_platform,
    runenv=runenv,
  }

  output {
    File output_file = run_checksex.output_file
  }
}
