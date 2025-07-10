version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/bwa/idx.wdl"

workflow bwa_build_idx {
    meta {
        author: "Eddie Belter"
        version: "0.1"
        description: "Build BWA index from FASTA reference"
    }

    input {
        String name
        File fasta # GZIP OK
        String docker
        Int cpu
        Int memory
    }

    RunEnv runenv = {
      "docker": docker,
      "cpu": cpu,
      "memory": memory,
      "disks": 20,
    }

    call idx.run_build_idx { input:
        name=name,
        fasta=fasta,
        runenv=runenv
    }

    output {
        File idx = run_build_idx.idx
    }
}
