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
  # Popular params:
  # --write-index
  command <<<
    samtools sort ~{params} -O ~{output_fmt} --write-index -o "~{output_bam}##idx##~{output_bam}.bai" ~{bam}
  >>>

  output {
    File output_bam = output_bam
    File output_bai = glob("*.bai")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
  }
}
