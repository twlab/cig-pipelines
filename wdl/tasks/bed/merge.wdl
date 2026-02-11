version development

import "../../structs/runenv.wdl"

task run_union_bedgraphs {
  meta {
    author: "Eddie Belter"
    version: "0.1"
    description: "Merge BedGraphs Summing the Coverage"
  }

  input {
    Array[File] bgs
    String output_bg
    RunEnv runenv
  }

  command <<<
    bedtools unionbedg  -i ~{sep=" " bgs} | awk 'OFS="\t" {sum=0; for (i=4; i<=NF; i++) sum+=$i; print $1,$2,$3,sum}' > ~{output_bg}
  >>>

  output {
    File bg = output_bg
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
  }
}
