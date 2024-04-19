version development

# From: https://github.com/vgteam/vg/wiki/Extracting-a-FASTA-from-a-Graph

import "../../structs/runenv.wdl"

task run_extract_ref {
  input {
     String name
     File gbz
     RunEnv runenv
  }

  String output_paths_list = "ref/~{name}.paths.txt"
  String output_fasta = "ref/~{name}.fasta"
  String output_fai = "ref/~{output_fasta}.fai"
  String output_dict = "ref/~{name}.dict"
  command <<<
    set -ex
    mkdir ref
    vg paths --list --xg ~{gbz} | grep -v _decoy | grep -v _random |  grep -v chrUn_ | grep -v chrEBV | grep -v chrM | grep ~{name} > ~{output_paths_list}
    vg paths --extract-fasta -x ~{gbz} --paths-file ~{output_paths_list} > ~{output_fasta}
    samtools faidx ~{output_fasta} -o ~{output_fai}
    samtools dict ~{output_fasta} -o ~{output_dict}
  >>>

  output {
    Directory path = "ref"
    File fasta = glob("~{output_fasta}")[0]
    File fai = glob("~{output_fasta}.fai")[0]
    File dict = glob("~{output_dict}")[0]
    File paths_list = glob("~{output_paths_list}")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "${runenv.memory} GB"
    #disks: runenv.disks
  }
}
