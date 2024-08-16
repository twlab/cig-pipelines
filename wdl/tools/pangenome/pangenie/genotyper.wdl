version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/misc/cat.wdl"
import "wdl/tasks/pangenie.wdl"
import "wdl/tasks/vcallers/utils.wdl"

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

  call utils.run_bgzip_and_index as bgzip_and_index_vcf { input:
    vcf=run_genotyper.vcf,
    runenv=runenv_pangenie_small,
  }

  output {
    File vcf = bgzip_and_index_vcf.vcf_gz
    File vcf_tbi = bgzip_and_index_vcf.vcf_tbi
    File histo = run_genotyper.histo
  }
}
