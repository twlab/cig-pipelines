version development

import "../../structs/runenv.wdl"
import "../../tasks/gatk/bqsr.wdl"

workflow gatk_bqsr {
    meta {
        author: "Eddie Belter"
        version: "0.1"
        description: "GATK BQSR"
    }

    input {
        File bam 
        Directory reference  # dir w/ fasta, fai, dict
        File known_sites     # vcf
        String docker = "broadinstitute/gatk:4.3.0.0"
        Int cpu = 4
        Int memory = 20
    }

    RunEnv runenv = {
      "docker": docker,
      "cpu": cpu,
      "memory": memory,
      "disks": 20,
    }

    call bqsr.run_bqsr { input:
        bam=bam,
        reference=reference,
        known_sites=known_sites,
        runenv=runenv
    }

    call bqsr.apply_bqsr { input:
        bam=bam,
        table=run_bqsr.table,
        reference=reference,
        runenv=runenv
    }

    output {
        File recal_bam = apply_bqsr.recal_bam
    }
}
