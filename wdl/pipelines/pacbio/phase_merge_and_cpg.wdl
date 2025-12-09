version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/pacbio/cpg.wdl"
import "wdl/tasks/pacbio/phase.wdl"
import "wdl/tasks/samtools/merge.wdl"
import "wdl/tasks/samtools/sort.wdl"

workflow phase_merge_call_cpg {
  meta {
    author: "Eddie Belter"
    version: "1.0"
  }

  input {
    Array[File] crams
    String reference
    File phased_reads_fof
    String pbcpg_params
    String output_cram_fn
    # Runtime (dockers, cpu, mem)
    String phasing_docker # has phase-reads script and samtools
    Int phasing_cpu
    Int phasing_memory
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

  call phase.run_phase_crams { input:
    crams=crams,
    reference=reference,
    phased_reads_fof=phased_reads_fof,
    runenv=phasing_runenv,
  }

  scatter(cram in run_phase_crams.phased_crams) {
    call sort.run_sort_cram { input:
      cram=cram,
      output_cram_fn=basename(cram, ".cram") + ".sorted.cram",
      reference=reference,
      runenv=sort_runenv,
    }
  }

  call merge.run_merge_crams { input:
    crams=run_sort_cram.sorted_cram,
    reference=reference,
    output_cram_fn=output_cram_fn,
    runenv=merge_runenv,
  }

  call cpg.run_cpg_cram { input:
    cram=run_merge_crams.merged_cram,
    crai=run_merge_crams.merged_crai,
    reference=reference,
    params=pbcpg_params,
    runenv=pbcpg_runenv,
  }
}
