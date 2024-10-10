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

task run_validate_vcf {
   input {
     String sample
     File vcf
     RunEnv runenv
   }

  command <<<
    bcftools view ~{vcf} > /dev/null
    RV=$?
    if [ $RV -eq 0 ]; then
      status="OK"
    else
      status="FAIL"
    fi 
    printf "%s\t%s\t%s" ~{sample} ~{basename(vcf)} ${status} | tee output
  >>>

  output {
    status
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks : runenv.disks
  }
}
