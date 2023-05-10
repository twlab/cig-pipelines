version development

import "../../structs/runenv.wdl"

task run_bwa_mem {
    input {
        String name
        Array[File] fastqs
        File idx
        RunEnv runenv
    }

    String bam = "${name}.bam"
    command <<<
        mkdir reference
        cd reference
        tar -xvf ~{idx}
        index_folder=$(ls)
        reference_fasta=$(ls | head -1)
        reference_folder=$(pwd)
        reference_index_path=$reference_folder/$reference_fasta
        cd ..

        bwa \
            mem \
            -t ~{runenv.cpu} \
            -K 320000000 \
            $reference_index_path \
            ~{fastqs[0]} \
            ~{default="" fastqs[1]} | \
            samtools view -hbS - > ~{bam}
    >>>

    runtime {
        docker: runenv.docker
        cpu : runenv.cpu
        memory : "${runenv.memory} GB"
    }

    output {
        File bam = "${bam}"
    }
}
