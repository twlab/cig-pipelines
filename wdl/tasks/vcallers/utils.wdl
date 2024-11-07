version development

import "../../structs/runenv.wdl"

# REQUIRES: bgzip and tabix
task run_bgzip_and_index {
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

task run_quick_validate_vcf {
  # this task use bcftools view to validate a vcf
  input {
    String input_sample
    File input_vcf
    Int exit_on_fail = 0
    RunEnv runenv
  }

  command <<<
    printf "Validate VCF: %s" ~{input_vcf}
    bcftools view ~{input_vcf} > /dev/null
    rv=$?
    test "${rv}" == "0" && status="PASS" || status="FAIL"
    printf "VCF status: %s" "${status}"
    printf "${status}" > status
    test ~{exit_on_fail} == 1 && exit "${rv}"
  >>>

  output {
    String sample = input_sample
    String vcf = input_vcf
    String status = read_string("status")
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks : runenv.disks
  }
}
