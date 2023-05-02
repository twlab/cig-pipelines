version development

# Split BAM into SAMS by CHRs

import "struct/runenv.wdl"
import "tools/samtools.wdl"

workflow samtools_stat {
    input {
        File bam
        String docker = "ebelter/samtools:1.15.1"
        Int cpu = 4
        Int memory = 24
        Int disks = 20
    }

    RunEnv runenv = {
      "docker": docker,
      "cpu": cpu,
      "memory": memory,
      "disks": disks,
    }

    call samtools.sort as sort { input:
        bam=bam,
        runenv=runenv,
    }

    call samtools.index as index { input:
        bam=sort.sort_bam,
        runenv=runenv,
    }

    call samtools.split as split { input:
        bam=sort.sorted_bam,
        bai=index.bai,
        runenv=runenv,
    }

    output {
        File sorted_bam = sort.sorted_bam
        File sorted_bai = index.bai
        String prefix = split.prefix
        File header = split.header
        Array[File] sams = split.sams
    }
}

# TODO move into tasks, generalize
#      move command into script, add to docker
task split{
    # expects pos sorted bam with index
    input {
        File bam
        File bai
        RunEnv runenv
    }

    Int samtools_cpu = runenv.cpu - 1
    Int samtools_mem = runenv.memory - 2

    String lbam = basename(bam)
    String lbai = "~{basename(bam)}.bai"

    String output_basename = sub(basename(bam), ".bam", "")
    command <<<
        echo LINK BAM AND BAI to execution directory
        ln -s ~{lbam} ~{bam} 
        ln -s ~{lbai} ~{bai} 

        set -ex
        echo AUTOSOMAL CHR $(date)
        for i in {1..22}
        do
            c="chr${i}"
            echo Chromosome: ${c}
            time samtools view ~{lbam} ${c} -M -O SAM -o "~{output_basename}.${c}.sam"
        done

        echo MINOR CHR $(date)
        samtools view -H ~{lbam} | grep '@SQ' | awk '{print $2}' | awk -F: '{print $2}' | sort -n > all_chr
        for i in {1..22}; do echo "chr${i}"; done | sort -n > autosomal_chr
        comm -23 all_chr autosomal_chr > minor_chr
        samtools view -H ~{lbam} | grep -f minor_chr | grep '@SQ' | sed 's/@SQ\t//' | while read sq; do sn=$(echo ${sq} | tr " " "\n" | grep SN | awk -F: '{print $2}'); ln=$(echo ${sq} | tr " " "\n" | grep LN | awk -F: '{print $2}'); echo -e "${sn}\t1\t${ln}"; done > minor_chr.bed
        time samtools view ~{lbam} -L minor_chr.bed -M -O SAM -o ~{output_basename}.minor.sam
        
        for f in *.sam
        do
            [[ -s ${f} ]] || rm -f ${f}
        done

        echo HEADER
        samtools view -H ~{lbam} -O SAM -o  ~{output_basename}.header
    >>>

    output {
        String prefix = output_basename
        File header = glob(output_basename+"*.header")[0]
        Array[File] sams = glob(output_basename+".*.sam")
    }

    runtime {
        cpu: runenv.cpu
        memory: "~{runenv.memory} GB"
        docker: runenv.docker
    }
}
