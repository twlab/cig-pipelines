version development

import "../../../structs/runenv.wdl"

task run_haplocheck {
  input {
    File vcf
    RunEnv runenv
  }

  # VCF & results dirs need to be absolute paths
  # filter vcf by chrM, putting file in path to pass to haplocheck
  command <<<
    set -x
    mkdir vcfs results
    results_dn=$(readlink -f results)
    vcfs_dn=$(readlink -f vcfs)
    vcf_bn=$(basename ~{vcf})
    tabix ~{vcf}
    bcftools view -r chrM -O z ~{vcf} > ${vcfs_dn}/${vcf_bn}
    /apps/haplocheck/cloudgene run haplocheck@1.3.2 --files ${vcfs_dn} --format VCF --output ${results_dn}  --threads ~{runenv.cpu}
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
