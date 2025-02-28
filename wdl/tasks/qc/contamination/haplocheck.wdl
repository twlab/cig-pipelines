version development

import "../../../structs/runenv.wdl"

task run_haplocheck {
  input {
    File bam
    File bai
    RunEnv runenv
  }

  # BAM & results dirs need to be absolute paths
  # BAI is required
  # Will need to update version in the future
  command <<<
    set -x
    mkdir data results
    data_dn=$(readlink -f data)
    results_dn=$(readlink -f results)
    ln ~{bam} ${data_dn}/
    ln ~{bai} ${data_dn}/
    time /apps/haplocheck/cloudgene run haplocheck@1.3.2 --files ${data_dn} --format bam --output ${results_dn}
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
