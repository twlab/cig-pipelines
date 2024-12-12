version development

import "../../structs/runenv.wdl"

task run_liftover {
  input {
    File paf
    File bed
    Int mapping_qual = 5
    Int alignment_length = 50000
    Int max_seq_divergence = 1
    RunEnv runenv
  }

  # Usage: paftools.js liftover [options] <aln.paf> <query.bed>
  # Options:
  #  -q INT    min mapping quality [5]
  #  -l INT    min alignment length [50000]
  #  -d FLOAT  max sequence divergence (>=1 to disable) [1]
  String bedfile = basename(paf, "paf") + "liftover.bed"
  command <<<
    paftools.js liftover -q ~{mapping_qual} -l ~{alignment_length} - d ~{max_seq_divergence} ~{paf} ~{bed} > ~{bedfile}
  >>>

  output {
    File bedfile = glob(bedfile)[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
    #disks: runenv.disks
  }
}
