version development

import "wdl/structs/runenv.wdl"

# https://github.com/vgteam/vg/wiki/Extracting-a-FASTA-from-a-Graph

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

  call run_paths_and_samtools { input:
    name=name,
    gbz=gbz,
    runenv=runenv,
  }

  output {
    File fasta = run_paths_and_samtools.fasta
    File fai = run_paths_and_samtools.fai
    File dict = run_paths_and_samtools.dict
  }
}

task run_paths_and_samtools {
  input {
     String name
     File gbz
     RunEnv runenv
  }

  String output_fasta = "~{name}.fasta"
  command <<<
    vg paths --extract-fasta -x ~{gbz} --paths-by ~{name} > ~{output_fasta}
    samtools dict ~{output_fasta} -o ~{name}.dict
    samtools faidx ~{output_fasta} -o ~{output_fasta}.fai
  >>>

  output {
    File fasta = "~{output_fasta}"
    File dict = "~{name}.dict"
    File fai = "~{output_fasta}.fai"
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "${runenv.memory} GB"
    #disks: runenv.disks
  }
}
