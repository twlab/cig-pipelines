version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/pacbio/cpg.wdl" as pbcpg

workflow pbcpg {
  meta {
    author: "Eddie Belter"
    version: "1.0"
  }

  input {
    File alignments       # BAM or CRAM
    File alignments_idx   # BAI or CRAI
    File reference_fasta  # for CRAM
    String pbcpg_params   # cpg params
    String pbcpg_docker
    Int pbcpg_cpu
    Int pbcpg_memory
  }

  RunEnv pbcpg_runenv = {
    "docker": pbcpg_docker,
    "cpu": pbcpg_cpu,
    "memory": pbcpg_memory,
    "disks": 20,
  }

  call pbcpg.run_cpg as cpg_calls { input:
    alignments=alignments,
    alignments_idx=alignments_idx,
    reference_fasta=reference_fasta,
    params=pbcpg_params,
    runenv=pbcpg_runenv,
  }

  output {
    File bed = cpg_calls.bed
    File bigwig = cpg_calls.bigwig
  }
}
