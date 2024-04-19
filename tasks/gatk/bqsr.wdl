version development

import "../../structs/runenv.wdl"

task run_bqsr {
    input {
        File bam
        Directory reference
        File known_sites
        RunEnv runenv
    }

    String bqsr_table = "recal_data.table"
    Int javamem = runenv.memory - 2
    command <<<
        reference_fasta=$(find "~{reference}" -name \*.fasta | head -1)
        reference_dict=$(find "~{reference}" -name \*.dict | head -1)
        gunzip -c ~{known_sites} > known_sites.vcf
        /gatk/gatk --java-options -Xmx4g IndexFeatureFile -I known_sites.vcf
        /gatk/gatk --java-options -Xmx~{javamem}g BaseRecalibrator \
           -I ~{bam} \
           -R $reference_fasta \
           --sequence-dictionary $reference_dict \
           --known-sites known_sites.vcf \
           -O ~{bqsr_table}
    >>>

    runtime {
        docker: runenv.docker
        cpu : runenv.cpu
        memory : "~{runenv.memory} GB"
    }

    output {
        File table = "~{bqsr_table}"
    }
}

task apply_bqsr {
    input {
        File bam
        Directory reference
        File table
        RunEnv runenv
    }

    String output_bam = sub(basename(bam), ".bam$", ".bqsr.bam")
    Int javamem = runenv.memory - 2
    command <<<
        reference_fasta=$(find ~{reference} -name \*.fasta | head -1)
        reference_dict=$(find ~{reference} -name \*.dict | head -1)
        /gatk/gatk --java-options -Xmx~{javamem}g ApplyBQSR \
            -I ~{bam} \
            -R $reference_fasta \
            --sequence-dictionary $reference_dict \
            --bqsr-recal-file ~{table} \
            -O ~{output_bam}
    >>>

    output {
        File recal_bam = "~{output_bam}"
    }

    runtime {
        docker: runenv.docker
        cpu : runenv.cpu
        memory : "~{runenv.memory} GB"
    }
}
