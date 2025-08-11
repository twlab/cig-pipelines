version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/pacbio/pbmm2.wdl"

workflow pbmm2_align {
  meta {
    author: "Eddie Belter"
    version: "1.0"
  }

  input {
    String sample
    File mmi
    File ref_fasta
    File pbmm2_input
    String pbmm2_params
    #Runtime
    String pbmm2_docker
    Int pbmm2_cpu
    Int pbmm2_memory
  }

  RunEnv pbmm2_runenv = {
    "docker": pbmm2_docker,
    "cpu": pbmm2_cpu,
    "memory": pbmm2_memory,
    "disks": 20,
  }

  call pbmm2.run_align_output_cram as align { input:
    sample=sample,
    bam=pbmm2_input,
    reference_mmi=mmi,
    reference_fasta=ref_fasta,
    params="~{pbmm2_params} --sample ~{sample}",
    runenv=pbmm2_runenv,
  }
}
