version development

import "../../structs/runenv.wdl"

task run_qualimap {
  input {
    File bam
    RunEnv runenv
  }

  # --outdir ?
  command <<<
    ln -s ~{bam} .
    qualimap bamqc -gd HUMAN -nt ~{runenv.cpu - 1} -outformat pdf -c --java-mem-size=~{runenv.memory - 2}G -bam ~{basename(bam)}
  >>>

  output {
    File results = glob("*_stats/genome_results.txt")[0]
    File report = glob("*_stats/report.pdf")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks : runenv.disks
  }
}
