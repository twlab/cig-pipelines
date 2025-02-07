version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/ancestory/snvstory.wdl"

workflow snvstory {
  meta {
    author: "Eddie Belter"
    version: "1.0"
  }

  input {
    File vcf
    Directory resource
    String genome_ver = "38"
    String mode = "WGS"
    String docker = "mgibio/snvstory:v1.1-buster"
    Int cpu = 2
    Int memory = 24
  }

  RunEnv snvstory_runenv = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
    "disks": 20,
  }

  call snvstory.run_igm_churchill_ancestry { input:
    input_vcf=vcf,
    resource=resource,
    genome_ver=genome_ver,
    mode=mode,
    runenv=snvstory_runenv,
  }
}
