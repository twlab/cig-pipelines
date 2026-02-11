version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/bed/genomecov.wdl"
import "wdl/tasks/bed/merge.wdl"
import "wdl/tasks/samtools/extract.wdl"

workflow phase_coverage {
  meta {
    author: "Eddie Belter"
    version: "1.0"
  }

  input {
    File haplotypes_tsv  # cram ref
    File sources_tsv     # cram ref
    String seqtools_docker
  }

  RunEnv seqtools_runenv = {
    "docker": seqtools_docker,
    "cpu": 1,
    "memory": 16,
    "disks": 20,
  }

  Array[Array[String]] haplotypes = read_tsv(haplotypes_tsv)
  # File names
  String hap1_bn = basename(haplotypes[0][0])
  String hap1_reads_fof = "~{hap1_bn}.reads.fof"
  String hap2_bn = basename(haplotypes[1][0])
  String hap2_reads_fof = "~{hap2_bn}.reads.fof"

  # Get haplotype reads from haplotype CRAM
   call extract.run_extract_reads as hap1_reads { input:
    sam=haplotypes[0][0],
    ref=haplotypes[0][1],
    output_fof=hap1_reads_fof,
    runenv=seqtools_runenv,
  }
  call extract.run_extract_reads as hap2_reads { input:
    sam=haplotypes[1][0],
    ref=haplotypes[1][1],
    output_fof=hap2_reads_fof,
    runenv=seqtools_runenv,
  }
 
  # Extract haplotype reads from source CRAM and run genomecov
  Array[Array[String]] crams_and_refs = read_tsv(sources_tsv)
  scatter ( source in crams_and_refs ) {
    # Extract haplotype reads from source CRAM and generate genomecov
     call genomecov.run_extract_and_genomecov as hap1_genomecov { input:
      sam=source[0],
      ref=source[1],
      reads_fof=hap1_reads.fof,
      runenv=seqtools_runenv,
    }
    call genomecov.run_extract_and_genomecov as hap2_genomecov { input:
      sam=source[0],
      ref=source[1],
      reads_fof=hap2_reads.fof,
      runenv=seqtools_runenv,
    }
  }

  # Merge BGs
  if ( length(crams_and_refs) > 1 ) {  # already "merged"
    call merge.run_union_and_sum_bedgraphs as hap1_merge_bgs { input:
      bgs=hap1_genomecov.bg,
      output_bg="~{hap1_bn}.bg",
      runenv=seqtools_runenv,
    }
    call merge.run_union_and_sum_bedgraphs as hap2_merge_bgs { input:
      bgs=hap2_genomecov.bg,
      output_bg="~{hap2_bn}.bg",
      runenv=seqtools_runenv,
    }
  } 

  output {
    File hap1_bg = select_first([hap1_merge_bgs.bg, hap1_genomecov.bg[0]])
    File hap2_bg = select_first([hap2_merge_bgs.bg, hap2_genomecov.bg[0]])
  }
}
