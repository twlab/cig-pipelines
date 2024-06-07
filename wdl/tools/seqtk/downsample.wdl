version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/seqtk/downsample.wdl"

workflow seqtk_downsample {
  input {
    Array[File] fastqs
    Float probability
    String docker
    Int cpu
    Int memory
    Int disks
  }

  RunEnv runenv = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
    "disks": disks,
  }

  call downsample.determine_fraction { input:
    fastq=fastqs[0],
    probability=probability,
    runenv=runenv,
  }

  scatter (i in range(length(fastqs))) {
    call downsample.run_downsample as ds { input:
      fastq=fastqs[i],
      fraction=determine_fraction.fraction,
      runenv=runenv,
    }
  }
  output {
    Array[File] output_fastqs = ds.output_fastq
  }
}
