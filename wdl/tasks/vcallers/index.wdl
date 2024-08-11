version development

import "../../structs/runenv.wdl"

# REQUIRES: bzip2 and tabix
task run_bzip2_and_index {
    input {
        File vcf
        RunEnv runenv
    }

  String vcf_bn = basename(vcf)
  String vcf_gz = "${vcf_bn}.gz"
  command <<<
    bgzip -c ~{vcf} > ~{vcf_gz}
    tabix -p vcf ~{vcf_gz}
  >>>

  output {
    File vcf_gz = glob("${vcf_gz}")[0]
    File vcf_tbi = glob("${vcf_gz}.tbi")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks : runenv.disks
  }
}
