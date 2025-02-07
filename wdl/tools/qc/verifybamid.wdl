version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/qc/verifybamid.wdl"

workflow verify_bam_id {
  input {
    File bam
    Directory index
    String sample
    Directory snvstory_resource
    String snvstory_genome_ver = "38"
    String snvstory_mode = "WGS"
    String docker
    Int cpu
    Int memory
    String snvstory_docker
    Int snvstory_cpu
    Int snvstory_memory
  }

  RunEnv runenv_verifybamid = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
    "disks": 20,
  }

  call verifybamid.run_verifybamid { input:
    input_bam=bam
    runenv=runenv_verifybamid,
  }
}
