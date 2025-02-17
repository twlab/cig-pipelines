version development

import "../../structs/runenv.wdl"

task run_qc_report {
  input {
    Array[File] ancestry_files
    Array[File] checksex_files
    Array[File] haplocheck_files
    Array[File] verifybamid_files
    String output_file = "qc_report.tsv"
    RunEnv runenv
  }

  command <<<
    smaht srbatch qc report
    python /apps/scripts/generate_qc_report.py \
      --ancestry-fof ~{write_lines(ancestry_files} \
      --sex-files ~{write_lines(checksex_files)} \
      --haplocheck-fof ~{write_lines(haplocheck_files} \
      --verifybamid-fof ~{write_lines(verifybamid_files} \
      --report-name ~{output_file}
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
