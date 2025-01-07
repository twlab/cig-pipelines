version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/samblaster.wdl"
import "wdl/tasks/bwa/idx.wdl"
import "wdl/tasks/picard/markdup.wdl"
import "wdl/tasks/samtools/index.wdl" as samtools_index
import "wdl/tasks/vcallers/deepvariant.wdl"

workflow markdup_and_deepvariant {
  input {
    File bam
    File idx         # tarred reference with FASTA, DICT, FAI
    String markduper # picard or samblaster
    String markdup_params
    String deepvariant_docker
    Int deepvariant_cpu
    Int deepvariant_memory
    String markdup_docker
    Int markdup_cpu
    Int markdup_memory
    String samtools_docker
    Int samtools_cpu
    Int samtools_memory
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

  RunEnv samtools_runenv = {
    "docker": samtools_docker,
    "cpu": samtools_cpu,
    "memory": samtools_memory,
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

  if ( markduper == "picard") {
    call markdup.run_markdup as picard_markdup { input:
      bam=bam,
      runenv=markduper_runenv,
    }
  }
  if ( markduper == "samblaster") {
    # Currently does not work as the BAM needs to be name sorted and in SAM format
    call samblaster.run_samblaster as samblaster_markdup { input:
      bam=bam,
      runenv=markduper_runenv,
    }
  } 

  File dedup_bam = select_first([picard_markdup.dedup_bam, samblaster_markdup.dedup_bam])

  call samtools_index.run_index { input:
    bam=dedup_bam,
    runenv=samtools_runenv,
  }

  call deepvariant.NEW_run_deepvariant as dv { input:
    sample=basename(dedup_bam, ".dedup.bam"),
    bam=dedup_bam,
    bai=run_index.bai,
    ref_fasta=reference.fasta,
    ref_fai=reference.fai,
    ref_dict=reference.dict,
    runenv=dv_runenv,
  }
}
