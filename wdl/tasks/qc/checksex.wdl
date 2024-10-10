version development

import "../structs/runenv.wdl"

task run_checksex {
  input {
    File input_vcf
    File input_vcf_tbi
    String seq_platform
    RunEnv runenv
  }

  String output_file = sub(basename(input_vcf), ".vcf.gz$", ".checksex.tsv")
  String input_vcf_bn = basename(input_vcf)
  String input_vcf_tbi_bn = basename(input_vcf_tbi)
  command <<<
    set -x
    ln -s ~{input_vcf} .
    ln -s ~{input_vcf_tbi} .
    smaht tools checksex ~{input_vcf_bn} ~{output_file} -p ~{seq_platform}
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
