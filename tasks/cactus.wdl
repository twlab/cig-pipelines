version development

import "../structs/runenv.wdl"

task minigraph { # create the minigraph
    input {
        File sequences
        String reference_name
        String out_name
        RunEnv runenv
    }

    String gfa = "${out_name + '.gfa.gz'}"
    # cactus-minigraph "${JOBSTORE}" "${SEQS}" "${GFA}" --reference "${REF}"
    command <<<
        cactus-minigraph \
            jobstore \
            ~{sequences} \
            ~{gfa} \
            ~{"--reference " + reference_name} \
            --binariesMode local \
            --maxMemory ~{runenv.memory - 2}G
    >>>

    output {
        File gfa = gfa
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: "~{runenv.memory} GB"
        disks: runenv.disks
    }
}

task graphmap { # map back to the minigraph
    input {
        File sequences
        File gfa
        String reference_name
        String out_name
        RunEnv runenv
    }

    String fasta = "${out_name + '.fa.gz'}"
    String paf = "${out_name + '.paf'}"
    String updated_sequences = basename(sequences)
    # cactus-graphmap "${JOBSTORE}" "${SEQS}" "${GFA}" "${PAF}" --outputFasta "${FASTA}" --reference "${REF}" --maxMemory 12G
    command <<<
        cp ~{sequences} ~{updated_sequences}
        cactus-graphmap \
            jobstore \
            ~{updated_sequences} \
            ~{gfa} \
            ~{paf} \
            ~{"--reference " + reference_name} \
            --outputFasta ~{fasta} \
            --binariesMode local \
            --maxMemory ~{runenv.memory - 2}G
    >>>

    output {
        File fasta = fasta
        File paf = paf
        File updated_sequences = updated_sequences
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: "~{runenv.memory} GB"
        disks: runenv.disks
    }
}

task graphmapsplit { # split the graph by chromosome
    input {
        File sequences
        File fasta
        File gfa
        File paf
        String reference_name
        String out_name
        RunEnv runenv
    }

    # cactus-graphmap-split "${JOBSTORE}" "${SEQS}" "${GFA}" "${PAF}" --outDir "${CHROMS}" --reference "${REF}" --maxMemory 12G
    command <<<
        echo ~{sequences} && \
        echo ~{fasta} && \
        sed -i 's#_MINIGRAPH_\t.*#_MINIGRAPH_\tfile://~{fasta}#' ~{sequences} && \
        cactus-graphmap-split \
            jobstore \
            ~{sequences} \
            ~{gfa} \
            ~{paf} \
            --outDir chroms \
            --reference ~{reference_name} \
            --binariesMode local \
            --maxMemory ~{runenv.memory - 2}G
    >>>

    output {
        Directory chroms = "chroms"
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: "~{runenv.memory} GB"
        disks: runenv.disks
    }
}

task align {
    input {
        Directory chroms
        File fasta
        String reference_name
        String out_name
        RunEnv runenv
    }

    String alignments = "chrom-alignments"
    # cactus-align "${JOBSTORE}" "${CHROMS}/chromfile.txt" "${ALIGNMENTS}" --batch --pangenome --reference "${REF}" --outVG 
    command <<<
        mv ~{chroms} ./
        find chroms/seqfiles -name \*.seqfile -exec sed -i 's#file:///.\+/execution/##' {} \;
        cactus-align \
            jobstore \
            chroms/chromfile.txt \
            ~{alignments} \
            --batch \
            ~{"--reference " + reference_name} \
            --pangenome \
            --outVG \
            --binariesMode local \
            --maxMemory ~{runenv.memory - 2}G
    >>>

    output {
        Directory alignments = alignments
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: "~{runenv.memory} GB"
        disks: runenv.disks
    }
}

task minigraphjoin {
    input {
        Directory alignments
        String reference_name
        String out_name
        RunEnv runenv
    }

    # cactus-graphmap-join "${JOBSTORE}" --vg "${ALIGNMENTS}/*.vg" --hal "${ALIGNMENTS}/*.hal" --outDir "${OUTDIR}" --outName "${NAME}" --reference "${REF}" --vcf --giraffe clip
    command <<<
        cactus-graphmap-join \
            jobstore \
            --vg ~{alignments}/*.vg \
            --hal ~{alignments}/*.hal \
            --outDir . \
            --outName ~{out_name} \
            ~{"--reference " + reference_name} \
            --vcf --giraffe \
            --binariesMode local \
            --maxMemory ~{runenv.memory - 2}G
    >>>

    output {
        # test.d2.dist test.d2.gbz test.d2.min
        # test.full.hal test.stats.tgz	
        # test.gfa.gz
        # test.raw.vcf.gz test.raw.vcf.gz.tbi
        # test.vcf.gz test.vcf.gz.tbi
        File dist = glob(out_name +"*dist")
        File gbz = glob(out_name +"*gbz")
        File gfa = glob(out_name +"*.gfa.gz")
        File hal = glob(out_name +"*hal")
        File min = glob(out_name +"*min")
        File stats = glob(out_name +"*stats.tgz")
        File raw_vcf = glob(out_name +".raw.vcf.gz")
        File raw_vcf_tbi = glob(out_name +".raw.vcf.gz.tbi")
        File vcf = glob(out_name +".vcf.gz")
        File vcf_tbi = glob(out_name +".vcf.gz.tbi")
    }

    runtime {
        docker: runenv.docker
        cpu: runenv.cpu
        memory: "~{runenv.memory} GB"
        disks: runenv.disks
    }
}
