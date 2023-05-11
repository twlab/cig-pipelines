version development

import "../../structs/runenv.wdl"

task run_haplotype_caller {
    input {
        File bam
        Directory reference
        RunEnv runenv
    }

    String output_vcf = basename(bam) + ".vcf"
    Int javamem = runenv.memory - 2
    command {
        reference_fasta=$(find ~{reference} -name \*.fasta | head -1)
        reference_dict=$(find ~{reference} -name \*.dict | head -1)
        samtools index ${bam}
        /gatk/gatk --java-options -Xmx~{javamem}g HaplotypeCaller \
           -I ~{bam} \
           -R $reference_fasta \
           --sequence-dictionary $reference_dict \
           -O ~{output_vcf}
    }

    runtime {
        docker: runenv.docker
        cpu : runenv.cpu
        memory : runenv.memory + " GB"
    }

    output {
        File vcf = "~{output_vcf}"
        File vcf_idx = "~{output_vcf}.idx"
    }
}
