version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/bwa/align.wdl"
import "wdl/tasks/bwa/idx.wdl"
import "wdl/tasks/deepvariant.wdl"
import "wdl/tasks/picard/markdup.wdl"
import "wdl/tasks/samtools.wdl"

workflow genome_wgs {
  meta {
      author: "Eddie Belter"
      version: "0.1"
      description: "Align and Call Variants Pipleine"
  }

  input {
      String sample
      Array[Array[File]] fastqs  # read fastqs
      File idx            # tarred BWA index
  }

  RunEnv runenv_idx = {
    "docker": "ebelter/linux-tk:latest",
    "cpu": 1,
    "memory": 4,
    "disks": 20,

  }
  call idx.run_untar_idx as reference { input:
      idx=idx,
      runenv=runenv_idx,
  }

  RunEnv runenv_bwa = {
    "docker": "ebelter/bwa:0.7.17",
    "cpu": 6,
    "memory": 36,
    "disks": 20,
  }

  scatter (i in range(length(fastqs))) {
      call align.run_bwa_mem as align { input:
          sample=sample,
          library=sample+"-lib"+i,
          fastqs=fastqs[i],
          reference=reference.path,
          runenv=runenv_bwa,
      }
  }

  RunEnv runenv_merge = {
    "docker": "ebelter/samtools:1.15.1",
    "cpu": 4,
    "memory": 20,
    "disks": 20,

  }

  call samtools.merge_bams as merge { input:
      sample=sample,
      bams=align.bam,
      runenv=runenv_merge,
  }

  RunEnv runenv_samtools = {
    "docker": "ebelter/samtools:1.15.1",
    "cpu": 1,
    "memory": 4,
    "disks": 20,
  }

  call samtools.sort as samtools_sort { input:
      bam=merge.merged_bam,
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

  call samtools.stat as stat { input:
      bam=markdup.dedup_bam,
      runenv=runenv_samtools,
  } 

  call samtools.index as index { input:
      bam=markdup.dedup_bam,
      runenv=runenv_samtools,
  } 

  RunEnv runenv = {
    "docker": "google/deepvariant:1.5.0", # "google/deepvariant:1.5.0-gpu"
    "cpu": 9,
    "memory": 48,
    "disks": 20,
  }

  call deepvariant.run_deepvariant as dv { input:
    sample=sample,
    bam=markdup.dedup_bam,
    bai=index.bai,
    reference=reference.path,
    runenv=runenv,
  }

  output {
      File bam = markdup.dedup_bam
      File bai = index.bai
      File stats = stat.stats
      File dedup_metrics = markdup.metrics
      File vcf = dv.vcf
  }
}
