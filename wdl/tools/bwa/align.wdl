version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/bwa/align.wdl"

workflow bwa_align {
    meta {
        author: "Eddie Belter"
        version: "0.1"
        description: "Align with BWA MEM"
    }

    input {
        String sample
        String library
        Array[File] fastqs
        Directory reference
        String docker = "mgibio/bwa:0.7.17"
        Int cpu = 8
        Int memory = 48
    }

    RunEnv runenv = {
      "docker": docker,
      "cpu": cpu,
      "memory": memory,
      "disks": 20,
    }

    call align.run_bwa_mem { input:
        sample=sample,
        library=library,
        fastqs=fastqs,
        reference=reference,
        runenv=runenv
    }

    output {
        File bam = run_bwa_mem.bam
    }
}
