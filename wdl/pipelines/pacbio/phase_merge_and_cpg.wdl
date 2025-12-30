version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/pacbio/cpg.wdl"
import "wdl/tasks/pacbio/phase.wdl"
import "wdl/tasks/samtools/merge.wdl" as samtools_merge
import "wdl/tasks/samtools/sort.wdl" as samtools_sort
import "wdl/tasks/samtools/stats.wdl" as samtools_stats

workflow phase_merge_call_cpg {
  meta {
    author: "Eddie Belter"
    version: "1.0"
  }

  input {
    Array[File] haplotype1_crams
    Array[File] haplotype2_crams
    String haplotype1_reference
    String haplotype2_reference
    String pbcpg_params
    String haplotype1_output_cram
    String haplotype2_output_cram
    # Runtime (dockers, cpu, mem)
    String phasing_docker # has phase-reads script and samtools
    Int phasing_cpu
    Int phasing_memory
    Int stats_cpu
    Int stats_memory
    Int sort_cpu
    Int sort_memory
    Int merge_cpu
    Int merge_memory
    String pbcpg_docker 
    Int pbcpg_cpu
    Int pbcpg_memory
  }

  RunEnv phasing_runenv = {
    "docker": phasing_docker,
    "cpu": phasing_cpu,
    "memory": phasing_memory,
    "disks": 20,
  }

  RunEnv stats_runenv = {
    "docker": phasing_docker,
    "cpu": stats_cpu,
    "memory": stats_memory,
    "disks": 20,
  }

  RunEnv sort_runenv = {
    "docker": phasing_docker,
    "cpu": sort_cpu,
    "memory": sort_memory,
    "disks": 20,
  }

  RunEnv merge_runenv = {
    "docker": phasing_docker,
    "cpu": merge_cpu,
    "memory": merge_memory,
    "disks": 20,
  }

  RunEnv pbcpg_runenv = {
    "docker": pbcpg_docker,
    "cpu": pbcpg_cpu,
    "memory": pbcpg_memory,
    "disks": 20,
  }

  scatter(haplotype_crams in zip(haplotype1_crams, haplotype2_crams)) {
    call phase.run_phase_crams { input:
      haplotype1_cram=haplotype_crams.left,
      haplotype2_cram=haplotype_crams.right,
      haplotype1_reference=haplotype1_reference,
      haplotype2_reference=haplotype2_reference,
      runenv=phasing_runenv,
    }

    call samtools_sort.run_sort_cram as haplotype1_sorted_cram { input:
      cram=run_phase_crams.haplotype1_phased_cram,
      output_cram_fn=basename(run_phase_crams.haplotype1_phased_cram, ".cram") + ".sorted.cram",
      reference=haplotype1_reference,
      runenv=sort_runenv,
    }

    call samtools_sort.run_sort_cram as hapltype2_sorted_cram { input:
      cram=run_phase_crams.haplotype2_phased_cram,
      output_cram_fn=basename(run_phase_crams.haplotype2_phased_cram, ".cram") + ".sorted.cram",
      reference=haplotype2_reference,
      runenv=sort_runenv,
    }
  }

  call samtools_merge.run_merge_crams as merge_haplotype1_crams { input:
    crams=haplotype1_sorted_cram.sorted_cram,
    reference=haplotype1_reference,
    output_cram_fn=haplotype1_output_cram,
    runenv=merge_runenv,
  }

  call samtools_merge.run_merge_crams as merge_haplotype2_crams { input:
    crams=hapltype2_sorted_cram.sorted_cram,
    reference=haplotype2_reference,
    output_cram_fn=haplotype2_output_cram,
    runenv=merge_runenv,
  }

  call samtools_stats.run_stats as haplotype1_samtools_stats { input:
    sam_file=merge_haplotype1_crams.merged_cram,
    reference=haplotype1_reference,
    runenv=stats_runenv,
  }

  call samtools_stats.run_stats as haplotype2_samtools_stats { input:
    sam_file=merge_haplotype2_crams.merged_cram,
    reference=haplotype2_reference,
    runenv=stats_runenv,
  }

  call cpg.run_cpg_cram as haplotype1_cpg { input:
    cram=merge_haplotype1_crams.merged_cram,
    crai=merge_haplotype1_crams.merged_crai,
    reference=haplotype1_reference,
    params=pbcpg_params,
    runenv=pbcpg_runenv,
  }

  call cpg.run_cpg_cram as haplotype2_cpg { input:
    cram=merge_haplotype2_crams.merged_cram,
    crai=merge_haplotype2_crams.merged_crai,
    reference=haplotype2_reference,
    params=pbcpg_params,
    runenv=pbcpg_runenv,
  }

  output {
    File haplotype1_cram = merge_haplotype1_crams.merged_cram
    File haplotype1_crai = merge_haplotype1_crams.merged_crai
    File haplotype1_stats = haplotype1_samtools_stats.stats
    File haplotype1_cpg_bed = haplotype1_cpg.bed
    File haplotype1_cpg_bigwig = haplotype1_cpg.bigwig
    File haplotype2_cram = merge_haplotype2_crams.merged_cram
    File haplotype2_crai = merge_haplotype2_crams.merged_crai
    File haplotype2_stats = haplotype2_samtools_stats.stats
    File haplotype2_cpg_bed = haplotype2_cpg.bed
    File haplotype2_cpg_bigwig = haplotype2_cpg.bigwig
  }
}
