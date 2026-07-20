version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/realign/abra2.wdl"

workflow abra2_realign {
  input {
    File bam
    File bai
    File targets
    File reference_fasta
    File reference_fai
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

  call abra2.run_realigner as realign { input:
    in_bam_file=bam,
    in_bam_index_file=bai,
    in_target_bed_file=targets,
    in_reference_file=reference_fasta,
    in_reference_index_file=reference_fai,
    runenv=runenv,
  }

  output {
    File realigned_bam = realign.indel_realigned_bam
    File realigned_bai = realign.indel_realigned_bam_index
  }
}
