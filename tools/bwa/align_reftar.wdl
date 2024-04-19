version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/bwa/align.wdl"
import "wdl/tasks/bwa/idx.wdl"

workflow bwa_align {
    meta {
        author: "Eddie Belter"
        version: "0.1"
        description: "Align with BWA MEM"
    }

    input {
        String name
        Array[File] fastqs
        File reference
        String docker = "ebelter/bwa:0.7.17"
        Int cpu = 8
        Int memory = 48
    }

    RunEnv runenv_untar_idx = {
      "docker": "ebelter/linux-tk:latest",
      "cpu": 1,
      "memory": 4,
      "disks": 20,
    }

    call idx.run_untar_idx as untarred_reference { input:
        idx=reference,
        runenv=runenv_untar_idx,
    }

    RunEnv runenv = {
      "docker": docker,
      "cpu": cpu,
      "memory": memory,
      "disks": 20,
    }

    call align.run_bwa_mem { input:
        name=name,
        fastqs=fastqs,
        reference=untarred_reference.path,
        runenv=runenv
    }

    output {
        File bam = run_bwa_mem.bam
        Directory reference_path = untarred_reference.path
    }
}
