version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/deepvariant.wdl"
import "wdl/tasks/samtools.wdl"
import "wdl/tasks/picard/markdup.wdl"
import "wdl/tasks/vg/giraffe.wdl"
import "wdl/tasks/vg/stats.wdl"
import "wdl/tasks/vg/surject.wdl"

workflow pangenome_wgs {
  meta {
    author: "Eddie Belter"
    version: "1.1"
    description: "Map WGS reads to the pangnome graph with giraffe, then surject alignments to call variants with run deepvariant."
  }

  input {
    String sample
    Array[File] fastqs
    File min
    File dist
    File gbz
    File reference_tar
    String docker = "quay.io/vgteam/vg:v1.48.0" #"quay.io/vgteam/vg@sha256:62a1177ab6feb76de6a19af7ad34352bea02cab8aa2996470d9d2b40b3190fe8"
    Int cpu
    Int memory
  }

  RunEnv runenv_giraffe = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
    "disks": 20,
  }

  call giraffe.run_giraffe { input:
    sample=sample,
    fastqs=fastqs,
    min=min,
    dist=dist,
    gbz=gbz,
    runenv=runenv_giraffe,
  }

  RunEnv runenv_vg = {
    "docker": docker,
    "cpu": 8,
    "memory": 64,
    "disks": 20,
  }

  call stats.run_stats as vg_stats { input:
    gam=run_giraffe.gam,
    runenv=runenv_vg,
  }

  call surject.run_surject { input:
    gam=run_giraffe.gam,
    sample=sample,
    library=sample+"-lib1",
    gbz=gbz,
    runenv=runenv_vg,
  }

  RunEnv runenv_samtools = {
    "docker": "ebelter/samtools:1.15.1",
    "cpu": 1,
    "memory": 4,
    "disks": 20,
  }

  call samtools.sort as samtools_sort { input:
      bam=run_surject.bam,
      runenv=runenv_samtools,
  } 

  RunEnv runenv_picard = {
    "docker": "ebelter/picard:2.27.4",
    "cpu": 4,
    "memory": 20,
    "disks": 20,
  }

  call markdup.run_markdup as markdup { input:
      sample=sample,
      bam=samtools_sort.sorted_bam,
      runenv=runenv_picard,
  }

  call samtools.stat as samtools_stat { input:
      bam=markdup.dedup_bam,
      runenv=runenv_samtools,
  } 

  call samtools.index as samtools_index { input:
      bam=markdup.dedup_bam,
      runenv=runenv_samtools,
  } 

  RunEnv runenv = {
    "docker": "google/deepvariant:1.6.0", # "google/deepvariant:1.6.0-gpu"
    "cpu": 9,
    "memory": 48,
    "disks": 20,
  }

  call deepvariant.run_deepvariant as dv { input:
    sample=sample,
    bam=markdup.dedup_bam,
    bai=samtools_index.bai,
    reference_tar=reference_tar,
    runenv=runenv,
  }

  output {
    File gam = run_giraffe.gam
    File gam_stats = vg_stats.stats
    File bam = markdup.dedup_bam
    File bai = samtools_index.bai
    File bam_stats = samtools_stat.stats
    File dv_vcf = dv.vcf
  }
}
