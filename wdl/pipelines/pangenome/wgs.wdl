version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/abra2.wdl"
import "wdl/tasks/bed/bedtools.wdl"
import "wdl/tasks/qc/fastqc.wdl"
import "wdl/tasks/freebayes.wdl"
import "wdl/tasks/gatk/realigner_target_creator.wdl"
import "wdl/tasks/pangenome/extract_ref.wdl"
import "wdl/tasks/picard/markdup.wdl"
import "wdl/tasks/samtools.wdl"
import "wdl/tasks/trimmers/fastp.wdl"
import "wdl/tasks/vcallers/deepvariant.wdl"
import "wdl/tasks/vg/giraffe.wdl"
import "wdl/tasks/vg/haplotypes.wdl"
import "wdl/tasks/vg/paths.wdl"
import "wdl/tasks/vg/stats.wdl"
import "wdl/tasks/vg/surject.wdl"

workflow pangenome_wgs {
  meta {
    author: "Eddie Belter"
    version: "1.2"
    description: "Map WGS reads to the pangnome graph with giraffe, then surject alignments bam, realign indels, and call variants with deepvariant."
  }

  input {
    String sample
    String reference_name
    Array[File] fastqs
    File gbz
    File hap
    Int targets_expansion_bases = 160
    Boolean generate_fastqc = false
    String? trimmer_name
    String? trimmer_params
    # dockers
    String abra2_docker 
    Int abra2_cpu 
    Int abra2_memory 
    String bedtools_docker
    Int bedtools_cpu
    Int bedtools_memory
    String deepvariant_docker
    Int deepvariant_cpu
    Int deepvariant_memory
    String? fastqc_docker
    Int? fastqc_cpu
    Int? fastqc_memory
    String freebayes_docker
    Int freebayes_cpu
    Int freebayes_memory
    String gatk_docker
    Int gatk_cpu
    Int gatk_memory
    String kmc_docker
    Int kmc_cpu
    Int kmc_memory
    String picard_docker
    Int picard_cpu
    Int picard_memory
    String samtools_docker
    Int samtools_cpu
    Int samtools_memory
    String? trimmer_docker
    Int? trimmer_cpu
    Int? trimmer_memory
    String vg_docker
    Int vg_cpu
    Int vg_memory
  }

  RunEnv abra2_renenv = {
    "docker": abra2_docker,
    "cpu": abra2_cpu,
    "memory": abra2_memory,
    "disks": 20,
  }

  RunEnv bedtools_runenv = {
    "docker": bedtools_docker,
    "cpu": bedtools_cpu,
    "memory": bedtools_memory,
    "disks": 20,
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

  RunEnv gatk_renenv = {
    "docker": gatk_docker,
    "cpu": gatk_cpu,
    "memory": gatk_memory,
    "disks": 20,
  }

  RunEnv giraffe_runenv = {
    "docker": vg_docker,
    "cpu": vg_cpu,
    "memory": vg_memory,
    "disks": 20,
  }

  RunEnv kmc_runenv = {
    "docker": kmc_docker,
    "cpu": kmc_cpu,
    "memory": kmc_memory,
    "disks": 20,
  }

  RunEnv picard_runenv = {
    "docker": picard_docker,
    "cpu": picard_cpu,
    "memory": picard_memory,
    "disks": 20,
  }

  RunEnv samtools_runenv = {
    "docker": samtools_docker,
    "cpu": samtools_cpu,
    "memory": samtools_memory,
    "disks": 20,
  }

  RunEnv vg_runenv = {
    "docker": vg_docker,
    "cpu": 8,
    "memory": 64,
    "disks": 20,
  }

  if ( generate_fastqc ) {
    RunEnv fastqc_runenv = {
      "docker": fastqc_docker,
      "cpu": fastqc_cpu,
      "memory": fastqc_memory,
      "disks": 20,
    }
    call fastqc.run_fastqc { input:
      seqfiles=fastqs,
      runenv=fastqc_runenv,
    }
  }

  if ( trimmer_name != "" ) {
    RunEnv trimmer_runenv = {
      "docker": trimmer_docker,
      "cpu": trimmer_cpu,
      "memory": trimmer_memory,
      "disks": 20,
    }
    if ( trimmer_name == "fastp") {
      call fastp.run_fastp as trimmer { input:
        fastqs=fastqs,
        params=trimmer_params,
        runenv=trimmer_runenv,
      }
    }
  }

  Array[File] trimmed_fastqs = select_first([trimmer.trimmed_fastqs, fastqs])

  call haplotypes.generate_kmers_with_kmc as kmc { input:
    sample=sample,
    fastqs=fastqs,
    runenv=kmc_runenv,
  }

  call giraffe.run_giraffe_haplotype_mode as run_giraffe { input:
    sample=sample,
    fastqs=fastqs,
    gbz=gbz,
    haplotypes=hap,
    kmers=kmc.kmers,
    runenv=giraffe_runenv,
  }

  call stats.run_stats as vg_stats { input:
    gam=run_giraffe.gam,
    runenv=vg_runenv,
  }

  call extract_ref.run_extract_ref as reference { input:
    name=reference_name,
    gbz=gbz,
    runenv=vg_runenv,
  }

  call surject.run_surject { input:
    gam=run_giraffe.gam,
    sample=sample,
    library=sample+"-lib1",
    gbz=gbz,
    paths_list=reference.paths_list,
    runenv=vg_runenv,
  }

  call samtools.sort as samtools_sort { input:
    bam=run_surject.bam,
    runenv=samtools_runenv,
  } 

  call freebayes.run_left_shift_bam as left_shift { input:
    in_bam_file=samtools_sort.sorted_bam,
    in_reference_file=reference.fasta,
    in_reference_index_file=reference.fai,
    runenv=freebayes_renenv,
  } 

  call samtools.index as left_shift_index { input:
    bam=left_shift.output_bam_file,
    runenv=samtools_runenv,
  } 

  call realigner_target_creator.run_realigner_target_creator as target_creator { input:
    in_bam_file=left_shift.output_bam_file,
    in_bam_index_file=left_shift_index.bai,
    in_reference_file=reference.fasta,
    in_reference_index_file=reference.fai,
    in_reference_dict_file=reference.dict,
    runenv=gatk_renenv,
  } 

  call bedtools.run_slop as expand_targets { input:
    bed_file=target_creator.output_target_bed_file,
    reference_fai=reference.fai,
    bases=targets_expansion_bases,
    runenv=bedtools_runenv,
  }

  call abra2.run_realigner as realign { input:
    in_bam_file=left_shift.output_bam_file,
    in_bam_index_file=left_shift_index.bai,
    in_target_bed_file=expand_targets.slopped_bed_file,
    in_reference_file=reference.fasta,
    in_reference_index_file=reference.fai,
    runenv=abra2_renenv,
  } 

  call markdup.run_markdup as picard_markdup { input:
    bam=realign.indel_realigned_bam,
    runenv=picard_runenv,
  }

  call samtools.index as samtools_index { input:
    bam=picard_markdup.dedup_bam,
    runenv=samtools_runenv,
  }

  call deepvariant.run_deepvariant as dv { input:
    sample=sample,
    bam=picard_markdup.dedup_bam,
    bai=samtools_index.bai,
    ref_fasta=reference.fasta,
    ref_fai=reference.fai,
    ref_dict=reference.dict,
    runenv=dv_runenv,
  }

  call samtools.stat as samtools_stat { input:
    bam=realign.indel_realigned_bam,
    runenv=samtools_runenv,
  } 

  output {
    File gam = run_giraffe.gam
    File gam_stats = vg_stats.stats
    File bam = picard_markdup.dedup_bam
    File bai = samtools_index.bai
    File bam_stats = samtools_stat.stats
    File vcf = dv.vcf
    File vcf_tbi = dv.vcf_tbi
    File gvcf = dv.gvcf
    File gvcf_tbi = dv.gvcf_tbi
  }
}
