version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/samtools/fastq.wdl"
import "wdl/tasks/pangenome/pangenie.wdl"

workflow pangenie_genotyper {
  input {
    File cram
    File reference
    File index
    String sample
    String params = ""
    String docker
    Int cpu
    Int memory
  }

  RunEnv runenv_cram2fastq = {
    "docker": docker,
    "cpu": 8,
    "memory": 32,
    "disks": 20,
  }

  RunEnv runenv_pangenie = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
    "disks": 20,
  }

  call fastq.run_sam_to_fastq { input:
    sam=cram,
    reference=reference,
    params="-F 0x900",
    runenv=runenv_cram2fastq,
  }

  call pangenie.run_genotyper { input:
    sample=sample,
    fastq=run_sam_to_fastq.fastq,
    index=index,
    params=params,
    runenv=runenv_pangenie,
  }

  output {
    File vcf = run_genotyper.vcf
    File vcf_tbi = run_genotyper.vcf_tbi
    File histo = run_genotyper.histo
  }
}
