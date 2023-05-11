version development

import "../../structs/runenv.wdl"

task run_build_idx {
    input {
        String name
        File fasta_gz
        RunEnv runenv
    }

    String fasta = "${name}.fasta"
    command <<<
        gunzip -c ~{fasta_gz} > ~{fasta}
        samtools dict ~{fasta} -o ~{fasta}.dict
        samtools faidx ~{fasta} -o ~{fasta}.fai
        bwa index -p ~{fasta} ~{fasta}
        tar cvvf ~{name}.tar ~{name + ".*"}
    >>>

    runtime {
        docker: runenv.docker
        cpu : runenv.cpu
        memory : "~{runenv.memory} GB"
    }

    output {
        File idx_tar = "${name}.tar"
    }
}
