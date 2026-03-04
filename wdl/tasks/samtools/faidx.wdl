version development

import "../../structs/runenv.wdl"

task extract_chromosome {
  input {
    File fasta
    File fai # or gzi
    String chr
    RunEnv runenv
  }

  String output_fn = "~{chr}.fasta"
  command <<<
    ln ~{fasta} ~{basename(fasta)}
    ln ~{fai} ~{basename(fai)}
    samtools faidx ~{fasta} ~{chr} -o ~{output_fn}
  >>>

  output {
    File chr_fasta = output_fn
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    #disks : select_first([runenv.disks,"local-disk 100 SSD"])
  }
}

task extract_chromosome_names {
  input {
    File fasta
    RunEnv runenv
  }

  String fai = "~{basename(fasta)}.fai"
  command <<<
    samtools faidx --fai-idx ~{fai} ~{fasta}
    cut -f1 ~{fai} | tee chromosomes
  >>>

  output {
    Array[String] names = read_lines("chromosomes")
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    #disks : select_first([runenv.disks,"local-disk 100 SSD"])
  }
}

task extract_chromosome_size {
  input {
    String chr
    File fai
    RunEnv runenv
  }

  command <<<
    awk '{if ($1 == "~{chr}"){print $2}}' ~{fai} | tee size
  >>>

  output {
    Int size = read_lines("size")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    #disks : select_first([runenv.disks,"local-disk 100 SSD"])
  }
}
