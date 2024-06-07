version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/samtools.wdl"
import "wdl/tasks/vg/map.wdl"
import "wdl/tasks/vg/stats.wdl"
import "wdl/tasks/vg/surject.wdl"

workflow pangenome_ataqseq {

  input {
    String sample
    Array[File] fastqs
    File gcsa
    File xg
    String docker = "quay.io/vgteam/vg:v1.48.0" #"quay.io/vgteam/vg@sha256:62a1177ab6feb76de6a19af7ad34352bea02cab8aa2996470d9d2b40b3190fe8"
    Int cpu = 12
    Int memory = 60
  }

  RunEnv runenv_map = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
    "disks": 20,
  }

  call map.run_map as map { input:
    sample=sample,
    fastqs=fastqs,
    gcsa=min,
    xg=xg,
    runenv=runenv_map,
  }

  RunEnv runenv_vg = {
    "docker": docker,
    "cpu": 8,
    "memory": 64,
    "disks": 20,
  }

  call stats.run_stats gam_stats { input:
    gam=map.gam,
    runenv=runenv_vg,
  }

  call surject.run_surject as surject { input:
    gam=map.gam,
    sample=sample,
    library=sample+"-lib1",
    gbz=gbz,
    runenv=runenv_vg,
  }

  RunEnv runenv_samtools = {
    "docker": docker, # vg 1.48.0 docker has samtools 1.10
    "cpu": 4,
    "memory": 20,
    "disks": 20,
  }

  call samtools.sort as samtools_sort { input:
    bam=surject.bam,
    runenv=runenv_samtools,
  }

  call samtools.index as samtools_index { input:
      bam=samtools_sort.sorted_bam,
      runenv=runenv_samtools,
  } 

  call samtools.stat as bam_stats { input:
    bam=samtools_sort.sorted_bam,
    runenv=runenv_samtools,
  }

  # TODO graph peak caller
  # call graph_peak_caller.run_grpah_peak_caller as peak_caller { input: }

  output {
    File gam = map.gam
    File gam_stats = gam_stats.stats
    File bam = samtools_sort.sorted_bam
    File bai = samtools_index.bai
    File bam_stats = bam_stat.stats
    #File peaks = peak_caller.peaks
  }
}
