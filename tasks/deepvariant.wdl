version development

import "../structs/runenv.wdl"

task run_deepvariant {
    input {
        String sample
        File bam
        File bai
        Directory reference
        RunEnv runenv
    }

    String output_vcf = "${sample}.vcf.gz"
    Int dv_cpu = runenv.cpu - 1
    command <<<
        ln ~{bam} ~{basename(bam)}
        ln ~{bai} ~{basename(bai)}
        reference_fasta=$(find ~{reference} -name \*.fasta)
        /opt/deepvariant/bin/run_deepvariant \
            --model_type=WGS \
            --ref=${reference_fasta} \
            --reads=~{basename(bam)} \
            --output_vcf=~{output_vcf} \
            --num_shards=~{dv_cpu}
    >>>

    output {
        String vcf = glob("${sample}.vcf.gz")[0]
        String vcf_tbi = glob("${sample}.vcf.gz.tbi")[0]
        String report = glob("${sample}*.html")[0]
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: runenv.memory + " GB"
        disks : runenv.disks
    }
}
