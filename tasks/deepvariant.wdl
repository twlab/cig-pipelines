version development

import "../structs/runenv.wdl"

task deep_variant {
    input {
        String name
        File bam
        File bai
        Directory reference
        RunEnv runenv
    }

    String output_vcf = "${name}.vcf.gz"
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
        String vcf = glob("${name}.vcf.gz")[0]
        String vcf_tbi = glob("${name}.vcf.gz.tbi")[0]
        String report = glob("${name}*.html")[0]
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: runenv.memory + " GB"
        #gpu_ount: 1
        disks : runenv.disks
    }
}
