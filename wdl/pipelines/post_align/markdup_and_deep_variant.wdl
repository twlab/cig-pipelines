version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/samblaster.wdl"
import "wdl/tasks/bwa/idx.wdl"
import "wdl/tasks/picard/markdup.wdl"
import "wdl/tasks/vcallers/deepvariant.wdl"

workflow markdup_and_deepvariant {
  input {
    String sample
    File bam
    File bai
    File idx         # tarred reference with FASTA, DICT, FAI
    String markduper # picard or samblaster
    String markdup_params
    String deepvariant_docker
    Int deepvariant_cpu
    Int deepvariant_memory
    String markdup_docker
    Int markdup_cpu
    Int markdup_memory
    String utils_docker
    Int utils_cpu
    Int utils_memory
  }

  RunEnv utils_runenv = {
    "docker": utils_docker,
    "cpu": utils_cpu,
    "memory": utils_memory,
    "disks": 20,
  }

  RunEnv markduper_runenv = {
    "docker": markdup_docker,
    "cpu": markdup_cpu,
    "memory": markdup_memory,
    "disks": 20,
  }

  RunEnv dv_runenv = {
    "docker": deepvariant_docker,
    "cpu": deepvariant_cpu,
    "memory": deepvariant_memory,
    "disks": 20,
  }

  call idx.run_untar_idx as reference { input:
    idx=idx,
    runenv=utils_runenv,
  }

  call markdup.run_markdup { input:
    bam=bam,
    runenv=markduper_runenv,
  }

  #call samtools.index as samtools_index { input:
  #  bam=markdup.decup_bam,
  #  runenv=samtools_runenv,
  #}

  call deepvariant.run_deepvariant as dv { input:
    sample=basename(run_markdup.dedup_bam, ".dedup.bam"),
    bam=run_markdup.dedup_bam,
    bai=bai,
    reference_path=reference.path,
    runenv=dv_runenv,
  }
}
