version development

import "../../structs/runenv.wdl"

task run_bwa_mem {
    input {
        String sample
        String library
        Array[File] fastqs
        Directory reference
        RunEnv runenv
    }

    String bam = "~{sample}.bam"
    Int bwa_cpu = runenv.cpu - 1
    command <<<
        reference_fasta=$(find ~{reference} -name \*.fasta)
        bwa \
            mem \
            -t ~{bwa_cpu} \
            -K 320000000 \
            -R '@RG\tID:~{library}\tLB:~{library}\tSM:~{sample}\tPL:illumina' \
            $reference_fasta \
            ~{fastqs[0]} \
            ~{default="" fastqs[1]} | \
            samtools view -hbS - > ~{bam}
    >>>

    runtime {
        docker: runenv.docker
        cpu : runenv.cpu
        memory : "~{runenv.memory} GB"
    }

    output {
        File bam = "~{bam}"
    }
}
