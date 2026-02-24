version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/pacbio/pbmm2.wdl"
import "wdl/tasks/pacbio/cpg.wdl"
import "wdl/tasks/samtools/stats.wdl" as samtools_stats

workflow align_and_cpg {
  meta {
    author: "Eddie Belter"
    version: "1.0"
  }

  input {
    String sample
    File pbmm2_input
    File reference_fasta
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

  RunEnv stats_runenv = {
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

  call pbmm2.run_align_output_cram as align { input:
    sample=sample,
    bam=pbmm2_input,
    reference_mmi=reference_mmi,
    reference_fasta=reference_fasta,
    params="~{pbmm2_params} --sample ~{sample}",
    runenv=pbmm2_runenv,
  }

  call samtools_stats.run_stats as samtools_stats { input:
    sam_file=align.aligned_cram,
    reference=reference_fasta,
    runenv=stats_runenv,
  }

  call cpg.run_cpg_cram as cpg_calls { input:
    cram=align.aligned_cram,
    crai=align.aligned_crai,
    reference=reference_fasta,
    params=pbcpg_params,
    runenv=pbcpg_runenv,
  }

  output {
    File cram = align.aligned_cram
    File crai = align.aligned_crai
    File stats = samtools_stats.stats
    File cpg_bed = cpg_calls.bed
    File cpg_bigwig = cpg_calls.bigwig
  }
}
