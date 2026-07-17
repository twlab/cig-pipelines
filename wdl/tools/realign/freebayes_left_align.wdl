version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/freebayes.wdl"

workflow freebayes_left_align {
  input {
    File bam
    File reference_fasta
    File reference_fai
    String docker
    Int cpu
    Int memory
  }

  RunEnv left_align_runenv = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
    "disks": 20,
  }

  call freebayes.run_left_shift_bam as left_align { input:
    in_bam_file=bam,
    in_reference_file=reference_fasta,
    in_reference_index_file=reference_fai,
    runenv=left_align_runenv,
  }

  output {
    File left_shifted_bam = left_align.output_bam_file
  }
}
