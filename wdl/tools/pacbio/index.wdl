version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/pacbio/index.wdl"

workflow pacbio_index {
  meta {
    author: "Eddie Belter"
    version: "0.1"
    description: "Build a Pacbio index from FASTA reference"
  }

  input {
    String mmi_name
    File fasta
    String params
    String pbmm2_docker
    Int pbmm2_cpu
    Int pbmm2_memory
  }

  RunEnv pbmm2_runenv = {
    "docker": pbmm2_docker,
    "cpu": pbmm2_cpu,
    "memory": pbmm2_memory,
    "disks": 20,
  }

  call index.run_index as build_index { input:
    mmi_name=mmi_name,
    fasta=fasta,
    params=params,
    runenv=pbmm2_runenv
  }
}
