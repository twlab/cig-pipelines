version development

import "../../structs/runenv.wdl"

task run_index {
  input {
    File fasta
    String mmi_name
    String params
    RunEnv runenv
  }

  # Alignment modes of --preset:
  #  SUBREAD     : -k 19 -w 10
  #  CCS or HiFi : -k 19 -w 10 -u
  #  ISOSEQ      : -k 15 -w 5  -u
  #  UNROLLED    : -k 15 -w 15
  command <<<
    pbmm2 index ~{params} -j ~{runenv.cpu} ~{fasta} ~{mmi_name}
  >>>

  output {
    File mmi = "~{mmi_name}"
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks : runenv.disks
  }
}
