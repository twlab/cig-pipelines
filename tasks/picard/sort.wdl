version development

import "../../structs/runenv.wdl"

task run_sort {
  input {
    File bam
    String sort_order
    RunEnv runenv
  }

  String output_bam = basename("~{bam}")
  command <<<
    java -Xmx~{runenv.memory - 1}g -jar /usr/picard/picard.jar SortSam \
      --INPUT ~{bam} \
      --OUTPUT ~{output_bam} \
      --SORT_ORDER ~{sort_order}
  >>>

  output {
    File output_bam = "~{output_bam}"
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks: runenv.disks
  }
}
