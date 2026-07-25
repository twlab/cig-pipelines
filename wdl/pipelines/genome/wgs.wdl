version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/bwa/align.wdl"
import "wdl/tasks/bwa/idx.wdl"
import "wdl/tasks/picard/markdup.wdl"
import "wdl/tasks/qc/fastqc.wdl"
import "wdl/tasks/gatk/realigner_target_creator.wdl"
import "wdl/tasks/realign/abra2.wdl"
import "wdl/tasks/realign/freebayes.wdl"
import "wdl/tasks/samtools.wdl"
import "wdl/tasks/trimmers/fastp.wdl"
import "wdl/tasks/vcallers/deepvariant.wdl"

struct StepsConf {
  Boolean realign_bam
}
struct Step {
  Boolean run
  String? params
  String docker
  Int cpu
  Int memory
}

workflow genome_wgs {
  meta {
      author: "Eddie Belter"
      version: "1.2"
      description: "Align with BWA, Left/Realign Indels then Call Variants with DeepVariant"
  }

  input {
    String sample
    Array[File] input_files
    File idx # tarred BWA index with DICT, FASTA, FAI
    String markdup_params = ""
    StepsConf steps_conf
    Step fastqc
    Step left_align_bam
    Step markdup_bam
    Int targets_expansion_bases = 160
    String? trimmer_name
    String? trimmer_params
    # resources
    String bwa_docker
    Int bwa_cpu
    Int bwa_memory
    String deepvariant_docker
    Int deepvariant_cpu
    Int deepvariant_memory
    String samtools_docker
    Int samtools_cpu
    Int samtools_memory
    String? trimmer_docker
    Int? trimmer_cpu
    Int? trimmer_memory
    String utils_docker
    Int utils_cpu
    Int utils_memory
    # realign resources
    String? abra2_docker
    Int? abra2_cpu
    Int? abra2_memory
    String? gatk_docker
    Int? gatk_cpu
    Int? gatk_memory
  }

  RunEnv bwa_runenv = {
    "docker": bwa_docker,
    "cpu": bwa_cpu,
    "memory": bwa_memory,
    "disks": 20,
  }

  RunEnv dv_runenv = {
    "docker": deepvariant_docker,
    "cpu": deepvariant_cpu,
    "memory": deepvariant_memory,
    "disks": 20,
  }

  RunEnv samtools_runenv = {
    "docker": samtools_docker,
    "cpu": samtools_cpu,
    "memory": samtools_memory,
    "disks": 20,
  }

  RunEnv utils_runenv = {
    "docker": utils_docker,
    "cpu": utils_cpu,
    "memory": utils_memory,
    "disks": 20,
  }

  call idx.run_untar_idx as reference { input:
    idx=idx,
    runenv=utils_runenv,
  }

  String input_ext = sub(input_files[0], "^.*\\.", "")

  if ( input_ext != "bam" ) {
    Array[String] fastqs = input_files
    if ( fastqc.run ) {
      RunEnv runenv_fastqc = {
        "docker": fastqc.docker,
        "cpu": fastqc.cpu,
        "memory": fastqc.memory,
        "disks": 20,
      }
      call fastqc.run_fastqc { input:
        seqfiles=fastqs,
        runenv=runenv_fastqc,
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

    call align.run_bwamem_with_sort as align { input:
      sample=sample,
      library=sample+"-lib1",
      rg_id=sample+"-lib1",
      platform_unit=sample+"-lib1",
      fastqs=trimmed_fastqs,
      idx_files=[reference.fasta, reference.amb, reference.ann, reference.bwt, reference.pac, reference.sa],
      runenv=bwa_runenv,
    }
  }

  File bam1 = select_first([align.bam, input_files[0]])

  if ( left_align_bam.run ) {
    RunEnv left_align_renenv = {
     "docker": left_align_bam.docker,
     "cpu": left_align_bam.cpu,
     "memory": left_align_bam.memory,
     "disks": 20,
    }

    call freebayes.run_left_align_bam as left_align { input:
      in_bam_file=bam1,
      in_reference_file=reference.fasta,
      in_reference_index_file=reference.fai,
      runenv=left_align_renenv,
    }
  }

  File bam2 = select_first([left_align.output_bam_file, bam1])

  if ( steps_conf.realign_bam ) {
    RunEnv abra2_renenv = {
      "docker": abra2_docker,
      "cpu": abra2_cpu,
      "memory": abra2_memory,
      "disks": 20,
    }
    RunEnv gatk_renenv = {
      "docker": gatk_docker,
      "cpu": gatk_cpu,
      "memory": gatk_memory,
      "disks": 20,
    }

    call samtools.index as bam2_idx { input:
      bam=bam2,
      runenv=samtools_runenv,
    }

    call realigner_target_creator.run_realigner_target_creator as target_creator { input:
      bam=bam2,
      bai=bam2_idx.bai,
      reference_fasta=reference.fasta,
      reference_fai=reference.fai,
      reference_dict=reference.dict,
      expand_bases=targets_expansion_bases,
      runenv=gatk_renenv,
    }

    call abra2.run_realigner as realign { input:
      in_bam_file=bam2,
      in_bam_index_file=bam2_idx.bai,
      in_reference_file=reference.fasta,
      in_reference_index_file=reference.fai,
      in_target_bed_file=target_creator.expanded_targets,
      runenv=abra2_renenv,
    }
  }

  File bam3 = select_first([realign.indel_realigned_bam, bam2])

  if ( markdup_bam.run ) {
    RunEnv markdup_bam_runenv = {
      "docker": markdup_bam.docker,
      "cpu": markdup_bam.cpu,
      "memory": markdup_bam.memory,
      "disks": 20,
    }

    call markdup.run_markdup as picard_markdup { input:
      bam=bam3,
      params=markdup_bam.params,
      runenv=markdup_bam_runenv,
    }
  }

  File bam4 = select_first([picard_markdup.dedup_bam, bam3])

  call samtools.index as bam4_idx { input:
    bam=bam4,
    runenv=samtools_runenv,
  }

  call samtools.stat as samtools_stat { input:
    bam=bam4,
    runenv=samtools_runenv,
  } 

  call deepvariant.run_deepvariant as dv { input:
    sample=sample,
    bam=bam4,
    bai=bam4_idx.bai,
    ref_fasta=reference.fasta,
    ref_fai=reference.fai,
    ref_dict=reference.dict,
    runenv=dv_runenv,
  }

  output {
    File bam = bam4
    File bai = bam4_idx.bai
    File bam_stats = samtools_stat.stats
    File vcf = dv.vcf
    File vcf_tvi = dv.vcf_tbi
  }
}
