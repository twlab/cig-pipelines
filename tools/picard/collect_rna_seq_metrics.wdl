version development

# Picard Metrics for RNA Seq

import "../../structs/runenv.wdl"
import "../../tasks/picard/collect_rna_seq_metrics.wdl"

workflow collect_rna_seq_metrics {
    input {
        File alignments     # [bam]
        File annotation     # [GTF GZ]
        String docker = "ebelter/picard:2.27.4"
        Int memory = 4 
        Int cpu = 20
    }

    RunEnv runenv = {
      "docker": docker,
      "cpu": cpu,
      "memory": memory,
      "disks": 20,
    }

    # Ribosomal Intervals
    call collect_rna_seq_metrics.create_ribosomal_intervals { input:
        alignments=alignments,
        annotation=annotation,
        runenv=runenv,
    }

    # REFLAT
    call collect_rna_seq_metrics.create_refflat { input:
        annotation=annotation,
        runenv=runenv,
    }

    # RNA Seq Metrics
    call collect_rna_seq_metrics.run_collect_rnaseq_metrics { input:
        alignments=alignments,
        ribosomal_intervals=create_ribosomal_intervals.intervals,
        refflat=create_refflat.refflat,
        runenv=runenv,
    }

    output {
        File metrics = run_collect_rnaseq_metrics.metrics
    }
}
