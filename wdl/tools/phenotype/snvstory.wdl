version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/snvstory.wdl"

workflow snvstory {
  meta {
    author: "Eddie Belter"
    version: "1.0"
  }

  input {
    Array[File] input_vcfs
    Directory resource
    String genome_ver = "38"
    String mode = "WGS"
    String docker = "mgibio/snvstory:v1.0-buster"
    Int cpu = 8
    Int memory = 64
  }

  RunEnv snvstory_runenv = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
    "disks": 20,
  }

  scatter(input_vcf in input_vcfs) {
    call snvstory.run_igm_churchill_ancestry { input:
      input_vcf=input_vcf,
      resource=resource,
      genome_ver=genome_ver,
      mode=mode,
      runenv=snvstory_runenv,
    }
  }
}
