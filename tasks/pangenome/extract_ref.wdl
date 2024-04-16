version development

# From: https://github.com/vgteam/vg/wiki/Extracting-a-FASTA-from-a-Graph

import "../wdl/structs/runenv.wdl"

task run_extract_ref {
  input {
     String name
     File gbz
     RunEnv runenv
  }

  String output_paths_list = "~{name}.paths.txt"
  String output_fasta = "~{name}.fasta"
  String output_fai = "~{output_fasta}.fai"
  String output_dict = "~{name}.dict"
  command <<<
    vg paths --list --xg ~{gbz} | grep -v _decoy | grep -v _random |  grep -v chrUn_ | grep -v chrEBV | grep -v chrM | grep ~{ref} > ~{output_paths_list}
    vg paths --extract-fasta -x ~{gbz} --paths-file ~{output_paths_listle} > ~{output_fasta}
    samtools faidx ~{output_fasta} -o ~{output_fai}
    samtools dict ~{output_fasta} -o ~{name}.dict
  >>>

  output {
    File fasta = glob("~{output_fasta}")[0]
    File fai = glob("~{output_fasta}.fai")[0]
    File dict = glob("~{output_dict}")[0]
    File paths_list = glob("~{paths_file}")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "${runenv.memory} GB"
    #disks: runenv.disks
  }
}
