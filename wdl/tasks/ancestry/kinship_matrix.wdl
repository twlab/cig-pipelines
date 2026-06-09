version development

import "../../structs/runenv.wdl"

task run_generate_kinship_matrix {
  input {
    Array[File] input_vcfs
    String output_file = "kinship"
    RunEnv runenv
  }

  command <<<
    set -x
    touch vcfs.fof
    for vcf in ~{sep=' ' input_vcfs}; do
      vcf_bn=$(basename "${vcf}")
      ln -s "${vcf}" "${vcf_bn}"
      echo "${vcf_bn}" >> vcfs.fof
    done
    python /apps/scripts/GenerateKinshipMatrix.py vcfs.fof ~{output_file}
  >>>

  output { # kinship.genome kinship.log kinship.nof kinship.nosex
    File output_file = glob("~{output_file}.genome")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "${runenv.memory} GB"
    #disks: runenv.disks
  }
}
