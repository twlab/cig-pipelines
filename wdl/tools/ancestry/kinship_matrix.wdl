version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/kinship_matrix.wdl"

workflow kinship_matrix {
  input {
    Array[File] vcfs
    String docker
    Int cpu
    Int memory
  }

  RunEnv runenv_km = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
    "disks": 20,
  }
  call kinship_matrix.run_generate_kinship_matrix { input:
    input_vcfs=vcfs,
    runenv=runenv_km,
  }

  output {
    File output_file = run_generate_kinship_matrix.output_file
  }
}
