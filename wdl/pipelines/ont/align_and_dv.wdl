version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/minimap2/align.wdl" as minimap2_align
import "wdl/tasks/samtools/sort.wdl" as samtools_sort
import "wdl/tasks/samtools/stats.wdl" as samtools_stats
import "wdl/tasks/vcallers/deepvariant.wdl"

workflow align_and_dv {
  meta {
    author: "Eddie Belter"
    version: "1.0"
  }

  input {
    String sample
    File reads
    File ref_fasta
    File ref_fai
    File ref_dict
    String minimap2_params
    String deepvariant_model_type = "ONT_R104"
    #Runtime
    String deepvariant_docker
    Int deepvariant_cpu
    Int deepvariant_memory
    String minimap2_docker
    Int minimap2_cpu
    Int minimap2_memory
    String samtools_docker
    Int samtools_cpu
    Int samtools_memory
    Int samtools_sort_cpu
    Int samtools_sort_memory
  }

  RunEnv dv_runenv = {
    "docker": deepvariant_docker,
    "cpu": deepvariant_cpu,
    "memory": deepvariant_memory,
    "disks": 20,
  }

  RunEnv minimap2_runenv = {
    "docker": minimap2_docker,
    "cpu": minimap2_cpu,
    "memory": minimap2_memory,
    "disks": 20,
  }

  RunEnv samtools_runenv = {
    "docker": samtools_docker,
    "cpu": samtools_cpu,
    "memory": samtools_memory,
    "disks": 20,
  }

  RunEnv samtools_sort_runenv = {
    "docker": samtools_docker,
    "cpu": samtools_sort_cpu,
    "memory": samtools_sort_memory,
    "disks": 20,
  }

  call minimap2_align.run_align as align { input:
    query=reads,
    target=ref_fasta,
    output_fn=sample + ".bam",
    params=minimap2_params,
    runenv=minimap2_runenv,
  }

  call samtools_sort.run_sort as sorted { input:
    bam=align.alignments,
    runenv=samtools_sort_runenv,
  } 

  call samtools_stats.run_stats { input:
    sam_file=sorted.sorted_bam,
    runenv=samtools_runenv,
  } 

  call deepvariant.run_deepvariant as dv { input:
    sample=sample,
    bam=sorted.sorted_bam,
    bai=sorted.sorted_bai,
    ref_fasta=ref_fasta,
    ref_fai=ref_fai,
    ref_dict=ref_dict,
    model_type=deepvariant_model_type,
    runenv=dv_runenv,
  }

  output {
    File final_bam = sorted.sorted_bam
    File final_bai = sorted.sorted_bai
    File stats = run_stats.stats
    File vcf = dv.vcf
    File vcf_tbi = dv.vcf_tbi
    File gvcf = dv.gvcf
    File gvcf_tbi = dv.gvcf_tbi
  }
}
