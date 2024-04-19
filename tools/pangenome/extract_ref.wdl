version development

import "wdl/tasks/pangenome/extract_:ref"
import "wdl/structs/runenv.wdl"

workflow pangenome_extract_ref {

  input {
    String name
    File gbz
    String docker = "quay.io/vgteam/vg:v1.48.0" #"quay.io/vgteam/vg@sha256:62a1177ab6feb76de6a19af7ad34352bea02cab8aa2996470d9d2b40b3190fe8"
    Int cpu = 4
    Int memory = 24
  }

  RunEnv runenv = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
    "disks": 20,
  }

  call extract_ref.run_extract_ref { input:
    name=name,
    gbz=gbz,
    runenv=runenv,
  }

  output {
    File fasta = run_extract_ref.fasta
    File fai = run_extract_ref.fai
    File dict = run_extract_ref.dict
    File paths_file = run_extract_ref.paths_file
  }
}
