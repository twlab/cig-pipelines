version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/metrics.wdl"

workflow seqlendist {
  input {
    String sample
    Array[File] seqfiles
    Array[String] labels
    Array[String] reports
    String bins = "lr"
  }

  RunEnv metrics_runenv = {
    "docker": "mgibio/rna-metrics:2.27.4",
    "cpu": 6,
    "memory": 24,
    "disks": 20,
  }

  call metrics.run_seqlendist { input:
    seqfiles=seqfiles,
    labels=labels,
    bins=bins,
    out=sample,
    reports=reports,
    runenv=metrics_runenv,
  }
}
