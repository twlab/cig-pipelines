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
        samtools dict ~{fasta} -o ~{name}.dict
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
        File idx = "${name}.tar"
    }
}

task run_untar_idx {
    input {
        File idx
        RunEnv runenv
		}

    command <<<
        mkdir ref
        cd ref
        tar -xvvfm ~{idx}
		>>>

    output {
        Directory path = "ref"
        File fasta = glob("ref/*.fasta")[0]
        File fai = glob("ref/*.fasta.fai")[0]
        File dict = glob("ref/*.dict")[0]
		}

    runtime {
        docker: runenv.docker
        cpu : runenv.cpu
        memory : runenv.memory +" GB"
    }
}
