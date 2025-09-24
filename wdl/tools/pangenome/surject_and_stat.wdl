version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/pangenome/extract_ref.wdl"
import "wdl/tasks/samtools.wdl"
import "wdl/tasks/vg/surject.wdl"

workflow surject_and_stat {
  meta {
    author: "Eddie Belter"
    version: "1.0"
    description: "Suject a GAM and run samtools stats on the resulting BAM."
  }

  input {
    String sample
    File gam
    String reference_name
    File gbz
    String samtools_docker
    Int samtools_cpu
    Int samtools_memory
    String vg_docker
    Int vg_cpu
    Int vg_memory
  }

  RunEnv samtools_runenv = {
    "docker": samtools_docker,
    "cpu": samtools_cpu,
    "memory": samtools_memory,
    "disks": 20,
  }

  RunEnv vg_runenv = {
    "docker": vg_docker,
    "cpu": 8,
    "memory": 64,
    "disks": 20,
  }

  call extract_ref.run_extract_ref as reference { input:
    name=reference_name,
    gbz=gbz,
    runenv=vg_runenv,
  }

  call surject.run_surject as vg_surject { input:
    gam=gam,
    sample=sample,
    library=sample+"-lib1",
    gbz=gbz,
    paths_list=reference.paths_list,
    runenv=vg_runenv,
  }

  call samtools.stat as samtools_stat { input:
    bam=vg_surject.bam,
    runenv=samtools_runenv,
  } 

  output {
    File bam = vg_surject.bam
    File stats = samtools_stat.stats
  }
}
