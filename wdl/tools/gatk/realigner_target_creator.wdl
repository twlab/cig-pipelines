version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/gatk/realigner_target_creator.wdl"

workflow gatk_realigner_target_creator {
  input {
    File bam
    File bai
    File reference_fasta
    File reference_fai
    File reference_dict
    Int expand_bases = 0
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

  call realigner_target_creator.run_realigner_target_creator as target_creator { input:
    bam=bam,
    bai=bai,
    reference_fasta=reference_fasta,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    expand_bases=expand_bases,
    runenv=runenv,
  }

  output {
    File intervals = target_creator.intervals
    File targets = target_creator.targets
    File expanded_targets = target_creator.expanded_targets
  }
}
