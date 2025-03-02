version development

import "../../../structs/runenv.wdl"

task run_haplocheck {
  input {
    String sample
    File bam
    File bai
    RunEnv runenv
  }

  # BAM & results dirs need to be absolute paths
  # BAI is required
  # Will need to update version in the future
  String bam_bn = basename(bam)
  command <<<
    set -x
    mkdir data results
    data_dn=$(readlink -f data)
    results_dn=$(readlink -f results)
    ln ~{bam} ${data_dn}/
    ln ~{bai} ${data_dn}/
    time /apps/haplocheck/cloudgene run haplocheck@1.3.2 --files ${data_dn} --format bam --output ${results_dn}
    sed -i 's/~{bam_bn}/~{sample}/' ${results_dn}/contamination/contamination.txt
    #sed -i 's/~{bam_bn}/~{sample}/' ${results_dn}/haplogroups/haplogroups.txt
  >>>

  output {
    File contamination_report = glob("results/contamination/contamination.txt")[0]
    File html_report = glob("results/report/report.html")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "${runenv.memory} GB"
    #disks: runenv.disks
  }
}
