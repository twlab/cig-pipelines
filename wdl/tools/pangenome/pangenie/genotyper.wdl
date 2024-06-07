version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/misc/cat.wdl"
import "wdl/tasks/pangenie.wdl"

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
    "docker": "ebelter/linux-tk:latest",
    "cpu": 1,
    "memory": 4,
    "disks": 10,
  }
  call cat.run_zcat as combined_fastq { input:
    files=fastqs,
    out="combined.fastq",
    runenv=runenv_cat,
  }

  RunEnv runenv_pangenie = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
    "disks": 20,
  }
  call pangenie.run_genotyper { input:
    fastq=combined_fastq.concatenated_file,
    index=index,
    sample=sample,
    runenv=runenv_pangenie,
  }

  output {
    File vcf = run_genotyper.vcf
    File histo = run_genotyper.histo
  }
}
