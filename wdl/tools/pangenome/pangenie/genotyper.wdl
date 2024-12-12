version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/misc/cat.wdl"
import "wdl/tasks/pangenome/pangenie.wdl"

workflow pangenie_genotyper {
  input {
    Array[File] fastqs
    Directory index
    String sample
    String docker
    Int cpu
    Int memory
  }

  RunEnv runenv_cat = {
    "docker": "mgibio/linux-tk:latest",
    "cpu": 1,
    "memory": 4,
    "disks": 10,
  }

  RunEnv runenv_pangenie = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
    "disks": 20,
  }

  RunEnv runenv_pangenie_small = {
    "docker": docker,
    "cpu": 1,
    "memory": 8,
    "disks": 20,
  }

  call cat.run_zcat as combined_fastq { input:
    files=fastqs,
    out="combined.fastq",
    runenv=runenv_cat,
  }

  call pangenie.run_genotyper { input:
    fastq=combined_fastq.concatenated_file,
    index=index,
    sample=sample,
    runenv=runenv_pangenie,
  }

  output {
    File vcf = run_genotyper.vcf_gz
    File vcf_tbi = run_genotyper.vcf_tbi
    File histo = run_genotyper.histo
  }
}
