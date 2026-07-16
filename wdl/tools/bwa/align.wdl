version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/bwa/align.wdl"

workflow bwa_align {
  meta {
    author: "Eddie Belter"
    version: "1.0"
    description: "Align FASTQS with BWA MEM to output a sorted BAM."
  }

  input {
    String sample
    Array[File] fastqs
    Array[File] idx_files
    String? library
    String? rg_id
    String? platform_unit
    String platform = "ILLUMINA"
    String docker
    Int cpu
    Int memory
  }

  String library_ = select_first([library, "~{sample}-lib1"])
  String rg_id_ = select_first([rg_id, library_])
  String platform_unit_ = select_first([platform_unit, rg_id_])

  RunEnv runenv = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
    "disks": 20,
  }

  call align.run_bwamem_with_sort as run_bwa { input:
    sample=sample,
    fastqs=fastqs,
    idx_files=idx_files,
    library=library_,
    rg_id=rg_id_,
    platform_unit=platform_unit_,
    platform=platform,
    runenv=runenv
  }

  output {
    File bam = run_bwa.bam
  }
}
