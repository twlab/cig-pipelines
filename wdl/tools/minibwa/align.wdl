version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/minibwa/align.wdl"

workflow minibwa_align {
  meta {
    author: "Eddie Belter"
    version: "1.0"
    description: "Align FASTQS with Mini-BWA to output a sorted BAM."
  }

  input {
    String sample
    Array[File] fastqs
    Array[File] idx_files # fasta l2b mbw
    String? library
    String? rg_id
    String? platform_unit
    String platform = "ILLUMINA"
    String params = ""
    Int bam_compression_level
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

  call align.run_minibwa {
    input:
      sample=sample,
      library=sample + "-lib1",
      rg_id=sample + "-lib1",
      platform_unit=sample + "-lib1",
      fastqs=fastqs,
      idx_files=idx_files,
      minibwa_params=params,
      bam_compression_level=bam_compression_level,
      runenv=runenv,
  }

  output {
    File bam = run_minibwa.bam
  }
}
