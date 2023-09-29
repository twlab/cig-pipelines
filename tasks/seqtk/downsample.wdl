version development

import "../../structs/runenv.wdl"

task run_downsample {
  input {
    File fastq
    Int fraction
    Int rng_seed = 37
    Boolean memory_saver = true
    RunEnv runenv
  }

  String output_fastq = basename("~{fastq}")
  command <<<
    seqtk sample ~{if (memory_saver) then "-2" else ""} -s~{rng_seed} ~{fastq} ~{fraction} | gzip - > ~{output_fastq}
  >>>

  output {
    File output_fastq = "~{output_fastq}"
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks: runenv.disks
  }
}

task determine_fraction {
  input {
    File fastq
    Float probability
    RunEnv runenv
  }

  command <<<
    zcat ~{fastq} | wc -l | awk '{ printf "%.0f\n", ($1/4) * ~{probability} }' > fraction
  >>>

  output {
    Int fraction = read_string("fraction")
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks: runenv.disks
  }
}
