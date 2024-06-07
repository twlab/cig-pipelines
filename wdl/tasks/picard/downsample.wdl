version development

import "../../structs/runenv.wdl"

task run_downsample {
  input {
    File bam
    String probability
    String params = ""
    RunEnv runenv
  }

  String output_bam = basename("~{bam}")
  String metrics_file = "~{output_bam}.metrics"
  command <<<
    java -Xmx~{runenv.memory - 1}g -jar /usr/picard/picard.jar DownsampleSam \
      --INPUT ~{bam} \
      --OUTPUT ~{output_bam} \
      --PROBABILITY ~{probability} \
      --METRICS_FILE ~{metrics_file} \
      ~{params}
  >>>

  output {
    File output_bam = "~{output_bam}"
    File metrics = "~{metrics_file}"
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks: runenv.disks
  }
}
