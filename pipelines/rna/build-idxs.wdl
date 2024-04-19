version 1.0

import "wdl/structs/runenv.wdl"

workflow build_idxs {
    input {
        File reference      # reference [fasta GZ]
        File annotation     # annotation [GTF GZ]
        String genome       # genome (e.g 'GRCh38')
        String anno_version # annotation version (e.g 'v24')
        String docker = "ebelter/rna-seq-pipeline:v1.2.11-nsp"
    }

    RunEnv runenv = {
      "docker": docker,
      "cpu": 4,
      "memory": 20,
      "disks": 10,
    }

    # Transcripts - build fasta and kallisto index
    call build_transcripts_fasta { input: # docker needs gffread
        reference=reference,
        annotation=annotation,
        runenv=runenv,
    }
    call build_transcripts_idx { input:
        reference=build_transcripts_fasta.transcripts,
        runenv=runenv,
    }

    # STAR
    RunEnv runenv_star = {
      "docker": docker,
      "cpu": 8,
      "memory": 64,
    }
    call build_star_idx { input:
        reference=reference,
        annotation=annotation,
        anno_version=anno_version,
        genome=genome,
        runenv=runenv_star,
    }

    # RSEM
    call build_rsem_idx { input:
        reference=reference,
        annotation=annotation,
        anno_version=anno_version,
        genome=genome,
        runenv=runenv,
    }

    output {
        File kallisto_index = build_transcripts_idx.index
        File rsem_index = build_rsem_idx.index
        File star_index = build_star_idx.index
        File transcripts = build_transcripts.transcripts
    }
}

task build_transcripts_fasta {
    input {
        File reference
        File annotation
        RunEnv runenv
    }

    command {
        $(which prep_transcripts.sh) \
            ~{reference} \
            ~{annotation}
    }

    output {
        File transcripts = "transcripts.fasta"
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: "~{select_first([runenv.memory,'8'])} GB"
    }
}

task build_transcripts_idx {
    input {
        File reference
        RunEnv runenv
    }

    command {
        $(which prep_kallisto.sh) ${reference}
    }

    output {
        File index = glob("*.idx")[0]
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: "~{select_first([runenv.memory,'8'])} GB"
    }
}

task build_star_idx {
    input {
        File reference
        File annotation
        String genome
        String anno_version
        RunEnv runenv
    }

    command {
        $(which prep_star.sh) \
            ~{reference} \
            ~{annotation} \
            ~{anno_version} \
            ~{genome} \
            ~{runenv.cpu}
    }

    output {
        File index = glob("*.tgz")[0]
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: "~{select_first([runenv.memory,'8'])} GB"
    }
}

task build_rsem_idx {
    input {
        File reference
        File annotation
        String genome
        String anno_version
        RunEnv runenv
    }

    command {
        $(which prep_rsem.sh) \
            ~{reference} \
            ~{annotation} \
            ~{anno_version} \
            ~{genome}
    }

    output {
        File index = glob("*.tgz")[0]
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: "~{select_first([runenv.memory,'8'])} GB"
    }
}
