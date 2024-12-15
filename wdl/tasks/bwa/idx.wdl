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
        tar -xvvf ~{idx}
        find . -type f -exec touch {} \;
		>>>

    output {
        Directory path = "ref"
        File dict = glob("ref/*.dict")[0]
        File fasta = glob("ref/*.fasta")[0]
        File fai = glob("ref/*.fasta.fai")[0]
        File sizes = glob("ref/*.sizes")[0]
        File amb = glob("ref/*.amb")[0]
        File ann = glob("ref/*.ann")[0]
        File bwt = glob("ref/*.bwt")[0]
        File pac = glob("ref/*.pac")[0]
        File sa = glob("ref/*.sa")[0]
		}

    runtime {
        docker: runenv.docker
        cpu : runenv.cpu
        memory : runenv.memory +" GB"
    }
}
