version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/picard/downsample.wdl"

workflow calculate_probability {
  input {
    File samtools_stat_file
    Int genome_size_gb
    Int coverage
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

  call downsample.run_calculate_probability { input:
    samtools_stat_file=samtools_stat_file,
    genome_size_gb=genome_size_gb,
    coverage=coverage,
    runenv=runenv,
  }
}
