version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/freebayes.wdl"
import "wdl/tasks/minimap2/align.wdl" as minimap2_align
import "wdl/tasks/samtools/index.wdl" as samtools_index
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
    #Runtime
    String deepvariant_docker
    Int deepvariant_cpu
    Int deepvariant_memory
    String freebayes_docker
    Int freebayes_cpu
    Int freebayes_memory
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

  RunEnv freebayes_renenv = {
    "docker": freebayes_docker,
    "cpu": freebayes_cpu,
    "memory": freebayes_memory,
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

  call freebayes.run_left_shift_bam as left_shift { input:
    in_bam_file=sorted.sorted_bam,
    in_reference_file=ref_fasta,
    in_reference_index_file=ref_fai,
    runenv=freebayes_renenv,
  }

  call samtools_index.run_index as left_shift_index { input:
    bam=left_shift.output_bam_file,
    runenv=samtools_runenv,
  } 

  call samtools_stats.run_stats { input:
    sam_file=left_shift.output_bam_file,
    runenv=samtools_runenv,
  } 

  call deepvariant.run_deepvariant as dv { input:
    sample=sample,
    bam=left_shift.output_bam_file,
    bai=left_shift_index.bai,
    ref_fasta=ref_fasta,
    ref_fai=ref_fai,
    ref_dict=ref_dict,
    runenv=dv_runenv,
  }

  output {
    File final_bam = left_shift.output_bam_file
    File final_bai = left_shift_index.bai
    File stats = run_stats.stats
    File vcf = dv.vcf
    File vcf_tbi = dv.vcf_tbi
    File gvcf = dv.gvcf
    File gvcf_tbi = dv.gvcf_tbi
  }
}
