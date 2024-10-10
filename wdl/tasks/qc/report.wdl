version development

import "../../structs/runenv.wdl"

task run_qc_report {
  input {
    Array[File] checksex_files
    Array[File] ancestory_files
    String output_file = "qc_report.tsv"
    RunEnv runenv
  }

  command <<<
    set -x
    python /apps/scripts/generate_qc_report.py --ancestry-files ~{sep=' ' ancestory_files} --sex-files ~{sep=' ' checksex_files} --report-name ~{output_file}
  >>>

  output {
    File output_file = glob("~{output_file}")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "${runenv.memory} GB"
    #disks: runenv.disks
  }
}
