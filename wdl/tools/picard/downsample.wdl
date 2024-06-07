version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/picard/downsample.wdl"
import "wdl/tasks/picard/helpers.wdl"
import "wdl/tasks/picard/sort.wdl"

workflow picard_downsample {
  input {
    String sample
    File bam
    Float probability
    String params
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

  call helpers.determine_sort_order as sort_order { input:
    bam=bam,
    runenv=runenv,
  }

  if ( ! sort_order.is_coordinate_sorted ) {
    call sort.run_sort as name_sort { input:
        bam=bam,
        sort_order="queryname",
        runenv=runenv,
    }
  }

  call downsample.run_downsample as ds { input:
    bam=select_first([name_sort.output_bam, bam]),
    probability=probability,
    params=params,
    runenv=runenv,
  }

  call sort.run_sort as coordinate_sort { input:
    bam=ds.output_bam,
    sort_order="coordinate",
    runenv=runenv,
  }

  output {
    File output_bam = coordinate_sort.output_bam
    File metrics = ds.metrics
  }
}
