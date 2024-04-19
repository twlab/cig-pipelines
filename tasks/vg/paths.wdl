version development

import "../../structs/runenv.wdl"

task run_paths_list_for_ref {
  input {
    File gbz
    String ref
    RunEnv runenv
  }

  # -L, --list   print (as a list of names, one per line) the path (or thread) names
  command <<<
    set -ex
    vg paths --list --xg ~{gbz} | grep -v _decoy | grep -v _random |  grep -v chrUn_ | grep -v chrEBV | grep -v chrM | grep ~{ref} > ~{ref}.list.txt
  >>>

  output {
    File paths_list = glob("~{ref}.list")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "${runenv.memory} GB"
    #disks: runenv.disks
  }
}
