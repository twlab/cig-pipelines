version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/bwa/align.wdl"
import "wdl/tasks/bwa/idx.wdl"
import "wdl/tasks/samtools/split.wdl"
import "wdl/tasks/wgsim.wdl"

workflow distortion_map {
    input {
      String sample
      File query_fasta
      File reference_idx      # tar file with ref fasta, fai, adn aligner index
      Array[String] chrs
      Float wgsim_base_error
      Int wgsim_out_distance
      Int wgsim_stdev
      Int wgsim_number_pairs
      Int wgsim_read1_length
      Int wgsim_read2_length
      Float wgsim_mutation_rate
      Float wgsim_fraction_indels
      Float wgsim_prob_indel_extentsion
      Int wgsim_seed
      # Resources
      String bwa_docker
      Int bwa_cpu
      Int bwa_memory
      String samtools_docker
      Int samtools_cpu
      Int samtools_memory
      String wgsim_docker
      Int wgsim_cpu
      Int wgsim_memory
      String utils_docker
      Int utils_cpu
      Int utils_memory
    }

  # RunEnvs in order of usage
  RunEnv utils_runenv = {
    "docker": utils_docker,
    "cpu": utils_cpu,
    "memory": utils_memory,
    "disks": 20,
  }

  RunEnv samtools_runenv = {
    "docker": samtools_docker,
    "cpu": samtools_cpu,
    "memory": samtools_memory,
    "disks": 20,
  }

  RunEnv wgsim_runenv = {
    "docker": wgsim_docker,
    "cpu": wgsim_cpu,
    "memory": wgsim_memory,
    "disks": 20,
  }

  RunEnv bwa_runenv = {
    "docker": bwa_docker,
    "cpu": bwa_cpu,
    "memory": bwa_memory,
    "disks": 20,
  }

  call idx.run_untar_idx as reference { input:
    idx=reference_idx,
    runenv=utils_runenv,
  }

  call split.run_split_by_chromosome as splitter { input:
    fasta=reference.fasta,
    fai=reference.fai,
    chrs=chrs,
    runenv=samtools_runenv,
  }

  scatter (chromosome_fasta in splitter.chromosome_fastas) {
    call wgsim.run_wgsim { input:
      fasta=chromosome_fasta,
      base_error=wgsim_base_error,
      out_distance=wgsim_out_distance,
      stdev=wgsim_stdev,
      number_pairs=wgsim_number_pairs,
      read1_length=wgsim_read1_length,
      read2_length=wgsim_read2_length,
      mutation_rate=wgsim_mutation_rate,
      fraction_indels=wgsim_fraction_indels,
      prob_indel_extentsion=wgsim_prob_indel_extentsion,
      seed=wgsim_seed,
      runenv=wgsim_runenv,
    }

    call align.run_bwa_mem as align_to_query { input:
      sample=sample,
      library=sample+"-lib1",
      fastqs=run_wgsim.simulated_fastqs,
      reference=reference.path,
      runenv=bwa_runenv,
    }
    call align.run_bwa_mem as align_to_ref { input:
      sample=sample,
      library=sample+"-lib1",
      fastqs=run_wgsim.simulated_fastqs,
      reference=reference.path,
      runenv=bwa_runenv,
    }
  }
}
