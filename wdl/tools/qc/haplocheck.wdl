version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/qc/contamination/haplocheck.wdl"

workflow haplocheck {
  input {
    String sample
    File bam
    File bai
    String docker
    Int cpu
    Int memory
  }

  RunEnv runenv = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
    "disks": 20,
  }

  call haplocheck.run_haplocheck { input:
    sample=sample,
    bam=bam,
    bai=bai,
    runenv=runenv,
  }

  output {
    File contamination_report = run_haplocheck.contamination_report
    File html_report = run_haplocheck.html_report
  }
}
