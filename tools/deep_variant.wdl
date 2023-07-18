version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/deepvariant.wdl"

workflow deep_variant {
    meta {
        author: "Eddie Belter"
        version: "0.1"
        description: "Call variants with Deep Variant"
    }

    input {
        String sample
        File bam
        File bai
        Directory reference
        String docker = "google/deepvariant:1.5.0" # "google/deepvariant:1.5.0-gpu"
        Int cpu = 9
        Int memory = 48
    }

    RunEnv runenv = {
      "docker": docker,
      "cpu": cpu,
      "memory": memory,
      "disks": 20,
    }

    call deepvariant.deep_variant { input:
        sample=sample,
        bam=bam,
        bai=bai,
        reference=reference,
        runenv=runenv
    }

    output {
        File vcf = deep_variant.vcf
        File vcf_tbi = deep_variant.vcf_tbi
        File repoort = deep_variant.report
    }
}
