version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/ancestory/snvstory.wdl"
import "wdl/tasks/misc/cat.wdl"
import "wdl/tasks/pangenome/pangenie.wdl"

workflow pangenie_genotyper_with_ancestory {
  input {
    String sample
    Array[File] fastqs
    File index
    String pangenie_params = ""
    File snvstory_resource
    String snvstory_genome_ver = "38"
    String snvstory_mode = "WGS"
    String docker
    Int cpu
    Int memory
    String snvstory_docker
    Int snvstory_cpu
    Int snvstory_memory
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

  RunEnv runenv_snvstory = {
    "docker": snvstory_docker,
    "cpu": snvstory_cpu,
    "memory": snvstory_memory,
    "disks": 20,
  }

  call cat.run_zcat as combined_fastq { input:
    files=fastqs,
    out="combined.fastq",
    runenv=runenv_cat,
  }

  call pangenie.run_genotyper { input:
    sample=sample,
    fastq=combined_fastq.concatenated_file,
    index=index,
    params=params,
    runenv=runenv_pangenie,
  }

  call snvstory.run_igm_churchill_ancestry { input:
    input_vcf=run_genotyper.vcf,
    resource=snvstory_resource,
    genome_ver=snvstory_genome_ver,
    mode=snvstory_mode,
    runenv=runenv_snvstory,
  }
}
