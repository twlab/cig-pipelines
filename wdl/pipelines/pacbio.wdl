version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/pacbio/pbmm2.wdl"
import "wdl/tasks/pacbio/cpg.wdl"
import "wdl/tasks/samtools/stats.wdl" as samtools_stats

workflow pacbio {
  meta {
    author: "Eddie Belter"
    version: "1.0"
  }

  input {
    String sample
    File pbmm2_input
    File reference_mmi
    String pbmm2_params
    String pbcpg_params
    #Runtime
    String pbmm2_docker
    Int pbmm2_cpu
    Int pbmm2_memory
    String pbcpg_docker 
    Int pbcpg_cpu
    Int pbcpg_memory
  }

  RunEnv pbmm2_runenv = {
    "docker": pbmm2_docker,
    "cpu": pbmm2_cpu,
    "memory": pbmm2_memory,
    "disks": 20,
  }

  RunEnv samtools_runenv = {
    "docker": pbmm2_docker,
    "cpu": 4,
    "memory": 16,
    "disks": 20,
  }

  RunEnv pbcpg_runenv = {
    "docker": pbcpg_docker,
    "cpu": pbcpg_cpu,
    "memory": pbcpg_memory,
    "disks": 20,
  }

  call pbmm2.run_align as align { input:
    sample=sample,
    bam=pbmm2_input,
    reference_mmi=reference_mmi,
    params="~{pbmm2_params} --sample ~{sample}",
    runenv=pbmm2_runenv,
  }

  call samtools_stats.run_stats as samtools_stats { input:
    sam_file=align.aligned_bam,
    runenv=samtools_runenv,
  }

  call cpg.run_cpg as cpg_scores { input:
    bam=align.aligned_bam,
    bai=align.aligned_bai,
    params=pbcpg_params,
    runenv=pbcpg_runenv,
  }
}
