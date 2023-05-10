version development

import "../../structs/runenv.wdl"
import "../../tasks/bwa/align.wdl"

workflow bwa_align {
    meta {
        author: "Eddie Belter"
        version: "0.1"
        description: "Align with BWA MEM"
    }

    input {
        String name
        Array[File] fastqs  # read fastqs
        File idx            # tarred BWA index
        String docker = "ebelter/bwa:0.7.17"
        Int cpu = 4
        Int memory = 20
    }

    RunEnv runenv = {
      "docker": docker,
      "cpu": cpu,
      "memory": memory,
    }

    call run_bwa_mem { input:
        name=name,
        fastqs=fastqs,
        idx=idx,
        runenv=runenv
    }

    output {
        File bam = run_bwa_mem.bam
    }
}
