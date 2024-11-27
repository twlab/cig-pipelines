version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/bwa/idx.wdl"
import "wdl/tasks/samtools/split.wdl"

workflow distortion_mapping {
    input {
        File query_fasta
        File reference_idx      # tar file with ref fasta, fai, adn aligner index
        String docker
        Int cpu
        Int memory
        String samtools_docker
        Int samtools_cpu
        Int samtools_memory
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

  call idx.run_untar_idx as reference { input:
    idx=idx,
    runenv=utils_runenv,
  }

  call split.run_split_by_chromosome as splitter { input:
    fasta=reference.fasta,
    fai=reference.fai,
    runenv=samtools_runenv,
  }

  #scatter (chromsome_fasta in splitter.chromsome_fastas) {
  #}
}
