version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/bed/genomecov.wdl"
import "wdl/tasks/pacbio/cov.wdl"
import "wdl/tasks/samtools/sort.wdl" as samtools_sort

workflow phase_coverage {
  meta {
    author: "Eddie Belter"
    version: "1.0"
  }

  input {
    File haplotype1_aligned_cram
    File haplotype1_aligned_reference
    File haplotype2_aligned_cram
    File haplotype2_aligned_reference
    File reference_aligned_cram
    File reference_aligned_ref
    String phasing_docker
    Int phasing_cpu
    Int phasing_memory
    String seqtools_docker
    Int seqtools_sort_cpu
    Int seqtools_sort_memory
    Int seqtools_genomecov_cpu
    Int seqtools_genomecov_memory
  }

  RunEnv phasing_runenv = {
    "docker": phasing_docker,
    "cpu": phasing_cpu,
    "memory": phasing_memory,
    "disks": 20,
  }

  RunEnv seqtools_sort_runenv = {
    "docker": seqtools_docker,
    "cpu": seqtools_sort_cpu,
    "memory": seqtools_sort_memory,
    "disks": 20,
  }

  RunEnv seqtools_genomecov_runenv = {
    "docker": seqtools_docker,
    "cpu": seqtools_genomecov_cpu,
    "memory": seqtools_genomecov_memory,
    "disks": 20,
  }

  call cov.extract_haplotype_reads_from_reference_aligned_cram as extract { input:
    haplotype1_aln_cram=haplotype1_aligned_cram,
    haplotype1_aln_ref=haplotype1_aligned_reference,
    haplotype2_aln_cram=haplotype1_aligned_cram,
    haplotype2_aln_ref=haplotype2_aligned_reference,
    reference_aln_cram=reference_aligned_cram,
    reference_aln_ref=reference_aligned_ref,
    runenv=phasing_runenv,
  }

  call samtools_sort.run_sort_cram as hap1_sort_cram { input:
    cram=extract.haplotype1_coverage_cram,
    output_cram_fn=basename(extract.haplotype1_coverage_cram, ".cram") + ".sorted.cram",
    reference=reference_aligned_ref,
    runenv=seqtools_sort_runenv,
  }

  call samtools_sort.run_sort_cram as hap2_sort_cram { input:
    cram=extract.haplotype2_coverage_cram,
    output_cram_fn=basename(extract.haplotype2_coverage_cram, ".cram") + ".sorted.cram",
    reference=reference_aligned_ref,
    runenv=seqtools_sort_runenv,
  }

  call genomecov.run_genomecov as hap1_genomecov { input:
    cram=hap1_sort_cram.sorted_cram,
    ref=reference_aligned_ref,
    runenv=seqtools_genomecov_runenv,
  }

  call genomecov.run_genomecov as hap2_genomecov { input:
    cram=hap2_sort_cram.sorted_cram,
    ref=reference_aligned_ref,
    runenv=seqtools_genomecov_runenv,
  }

  output {
    File haplotype1_cram = hap1_sort_cram.sorted_cram
    #File haplotype2_cram = 
    File haplotype1_cov_bed = hap1_genomecov.bed
    #File haplotype2_cov_bed = 
    File phase_coverage_stats = extract.stats
  }
}
