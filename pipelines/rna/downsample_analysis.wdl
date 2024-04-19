version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/misc/gunzip.wdl"
import "wdl/tasks/minimap2/lrna.wdl"
import "wdl/tasks/picard/collect_rna_seq_metrics.wdl"
import "wdl/tasks/seqtk/downsample.wdl"
import "wdl/tasks/rna/stringtie.wdl"

workflow downsample_analysis {
  input {
    String sample
    File fastq
    String read_type
    String strandedness
    Array[Float] probabilites
    Directory reference
		File annotation
		File junctions
  }

  RunEnv mm2_runenv = {
    "docker": "mgibio/minimap2:2.26",
    "cpu": 8,
    "memory": 48,
    "disks": 20,
  }

  RunEnv rnametrics_runenv = {
    "docker": "mgibio/rna-metrics:2.27.4",
    "cpu": 4,
    "memory": 20,
    "disks": 20,
  }

  RunEnv seqtk_runenv = {
    "docker": "mgibio/seqtools:latest",
    "cpu": 4,
    "memory": 16,
    "disks": 20,
  }

  RunEnv stringtie_runenv = {
    "docker": "mgibio/stringtie2:v2.0-focal",
    "cpu": 6,
    "memory": 24,
    "disks": 20,
  }

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

  # Downsample, align, collect rnaseq metrics, stringtie
  scatter (probability in probabilites) {
    call downsample.run_downsample { input:
      fastq=fastq,
      fraction=probability,
      runenv=seqtk_runenv,
    } # output_fastq

    call lrna.run_minimap2 as minimap2_align { input:
      fastq=run_downsample.output_fastq,
      reference=reference,
      junctions=junctions,
      output_prefix=sample,
      read_type=read_type,
      runenv=mm2_runenv,
    } # output_bam - primary only, coord sorted

    call stringtie.run_stringtie_denovo as stringtie_denovo { input:
      sample=sample,
      bam=minimap2_align.output_bam,
      strandedness=strandedness,
      params="-L -c 1.5 -f 0.05",
      runenv=stringtie_runenv,
    }

    call stringtie.run_stringtie_quantification as stringtie_quant { input:
      sample=sample,
      bam=minimap2_align.output_bam,
      annotation=annotation_gunzip.uncompressed,
      strandedness=strandedness,
      params="-L -c 1.5 -f 0.05",
      runenv=stringtie_runenv,
    }
  }
}
