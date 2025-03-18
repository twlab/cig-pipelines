version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/qc/contamination/haplocheck.wdl"

workflow haplocheck {
  input {
    File bams_tsv
    String haplocheck_docker
    Int haplocheck_cpu
    Int haplocheck_memory
  }

  RunEnv runenv_haplocheck = {
    "docker": haplocheck_docker,
    "cpu": haplocheck_cpu,
    "memory": haplocheck_memory,
    "disks": 20,
  }

  Array[Array[String]] bams = read_tsv(bams_tsv)
  scatter(info in bams) {
    String rg_id = info[0]
    String bam = info[1]
    String bai = info[2]

    call haplocheck.run_haplocheck { input:
      sample=rg_id,
      bam=bam,
      bai=bai,
      runenv=runenv_haplocheck,
    }
  }
}
