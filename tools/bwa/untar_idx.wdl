version development

import "../../structs/runenv.wdl"
import "../../tasks/bwa/idx.wdl"

workflow untar_idx {
    meta {
        author: "Eddie Belter"
        version: "0.1"
        description: "Untar a BWA Index"
    }

    input {
        File idx
        String docker = "ebelter/linux-tk:latest"
        Int cpu = 1
        Int memory = 4
    }

    RunEnv runenv = {
      "docker": docker,
      "cpu": cpu,
      "memory": memory,
      "disks": 20,
    }

    call idx.run_untar_idx as untarred_idx { input:
        idx=idx,
        runenv=runenv
    }

    output {
        Directory path = untarred_idx.path
        File fasta = untarred_idx.fasta
        File fai = untarred_idx.fai
        File dict = untarred_idx.dict
		}
}
