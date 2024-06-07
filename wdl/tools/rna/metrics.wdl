version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/metrics.wdl"
import "wdl/tasks/picard/collect_rna_seq_metrics.wdl"
import "wdl/tasks/picard/collect_insert_size_metrics.wdl"

workflow rna_metrics {
  input {
    String sample
    File fastq
    Array[String] seqlendist_reports
    String seqlendist_bin = "lr"
    File bam
    Array[String] lendist_reports
    File annotation
  }

  RunEnv metrics_runenv = {
    "docker": "mgibio/rna-metrics:2.27.4",
    "cpu": 6,
    "memory": 24,
    "disks": 20,
  }

  call metrics.run_seqlendist { input:
    seqfiles=[fastq],
    bin=seqlendist_bin,
    out=sample,
    reports=seqlendist_reports,
    runenv=metrics_runenv,
  }

  call collect_rna_seq_metrics.create_ribosomal_intervals { input:
    alignments=bam,
    annotation=annotation,
    runenv=metrics_runenv,
  }

  call collect_rna_seq_metrics.create_refflat { input:
    annotation=annotation,
    runenv=metrics_runenv,
  }

  call collect_rna_seq_metrics.run_collect_rnaseq_metrics { input:
    alignments=bam,
    ribosomal_intervals=create_ribosomal_intervals.intervals,
    refflat=create_refflat.refflat,
    runenv=metrics_runenv,
  }

  call collect_insert_size_metrics.run_collect_insert_size_metrics { input:
    alignments=bam,
    runenv=picard_runenv,
    runenv=metrics_runenv,
  }
}
