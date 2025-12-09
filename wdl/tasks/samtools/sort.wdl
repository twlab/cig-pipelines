version development

import "../../structs/runenv.wdl"

task run_sort {
  input {
    File bam
    String output_fmt = "BAM"
    String params = "" 
    RunEnv runenv
  }

  String output_bam = basename(bam)
  command <<<
    mkdir tmpsort
    trap rm -rf tmpsort EXIT
    samtools sort ~{params} -@ ~{runenv.cpu} -O ~{output_fmt} --write-index -o "~{output_bam}##idx##~{output_bam}.bai" ~{bam} -T tmpsort/sorted.nnnn.bam
  >>>

  output {
    File sorted_bam = output_bam
    File sorted_bai = glob("*.bai")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
  }
}

task run_sort_cram {
  input {
    File cram
    File reference
    String params = "" # --write-index
    String output_cram_fn
    RunEnv runenv
  }

  command <<<
    set -e
    mkdir tmpsort
    trap 'rm -rf tmpsort' EXIT
    samtools sort ~{params} --reference ~{reference} -@ ~{runenv.cpu} -O CRAM -o ~{output_cram_fn} ~{cram} -T tmpsort/sorted.nnnn.cram
  >>>

  output {
    File sorted_cram = "~{output_cram_fn}"
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
  }
}
