version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/bwa/idx.wdl"
import "wdl/tasks/samtools/split.wdl"
import "wdl/tasks/misc/tar.wdl"

workflow bwa_build_idx_by_chr {
  meta {
    author: "Eddie Belter"
    version: "0.1"
    description: "Build BWA index from FASTA reference"
  }

  input {
    String name
    File ref_idx
    String bwa_docker
    Int bwa_cpu
    Int bwa_memory
  }

  RunEnv bwa_runenv = {
    "docker": bwa_docker,
    "cpu": bwa_cpu,
    "memory": bwa_memory,
    "disks": 20,
  }

  # Untar the BWA Index
  call idx.run_untar_idx as ref { input:
    idx=ref_idx,
    runenv=bwa_runenv,
  }

  # Split the FASTA by Chromosome
  call split.run_split_by_chromosome as splitter { input:
    fasta=ref.fasta,
    fai=ref.fai,
    chrs=ref.chromosomes, #select_first([chrs, ref.chromosomes]),
    runenv=bwa_runenv,
  }

  # Build a BWA Index for Each Chromosome
  scatter (chromosome_fasta in splitter.chromosome_fastas) {
    call idx.run_build_idx { input:
      name=basename(chromosome_fasta, ".fasta"),
      fasta_gz=chromosome_fasta,
      runenv=bwa_runenv
    }
  }

  # TAR the CHR IDXs
  call tar.run_tar as tarred_indexes { input:
    name="~{name}.chr-bwa-idx",
    files=run_build_idx.idx,
    runenv=bwa_runenv,
  }

  output {
    File idx = tarred_indexes.tar_file
  }
}
