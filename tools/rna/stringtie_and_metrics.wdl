version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/misc/gunzip.wdl"
import "wdl/tasks/rna/stringtie.wdl"
import "wdl/tasks/picard/collect_rna_seq_metrics.wdl"
import "wdl/tasks/picard/collect_insert_size_metrics.wdl"

workflow stringtie_and_metrics {
  input {
    String sample
    File bam
    Boolean is_paired = false
    File annotation
    String strandedness # forward reverse unstranded
    Boolean denovo_mode = true
    String denovo_params = ""
    Boolean quant_mode = true
    String quant_params = ""
  }

  RunEnv stringtie_runenv = {
    "docker": "mgibio/stringtie2:v2.0-focal",
    "cpu": 6,
    "memory": 24,
    "disks": 20,
  }

  if ( denovo_mode ) {
    call stringtie.run_stringtie_denovo { input:
      sample=sample,
      bam=bam,
      strandedness=strandedness,
      params=denovo_params,
      runenv=stringtie_runenv,
    }
  }

  if ( quant_mode ) {
    RunEnv linuxtk_runenv = {
      "docker": "ebelter/linux-tk:latest",
      "cpu": 1,
      "memory": 4,
      "disks": 20,
    }
    call gunzip.run_gunzip as annotation_gunzip { input:
      compressed=annotation,
      runenv=linuxtk_runenv,
    }

    call stringtie.run_stringtie_quantification { input:
      sample=sample,
      bam=bam,
      annotation=annotation_gunzip.uncompressed,
      strandedness=strandedness,
      params=quant_params,
      runenv=stringtie_runenv,
    }
  }

  RunEnv picard_runenv = {
    "docker": "ebelter/picard:2.27.4",
    "cpu": 6,
    "memory": 24,
    "disks": 20,
  }

  call collect_rna_seq_metrics.create_ribosomal_intervals { input:
    alignments=bam,
    annotation=annotation,
    runenv=picard_runenv,
  }

  call collect_rna_seq_metrics.create_refflat { input:
    annotation=annotation,
    runenv=picard_runenv,
  }

  call collect_rna_seq_metrics.run_collect_rnaseq_metrics { input:
    alignments=bam,
    ribosomal_intervals=create_ribosomal_intervals.intervals,
    refflat=create_refflat.refflat,
    runenv=picard_runenv,
  }

  if ( is_paired ) {
    call collect_insert_size_metrics.run_collect_insert_size_metrics { input:
        alignments=bam,
        runenv=picard_runenv,
    }
  }
}
