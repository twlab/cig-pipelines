version development

import "../../structs/runenv.wdl"

task run_build_idx {
  # Create an BWA index TAR with support files.
  # CHM13v2.0.dict         samtools dict
  # CHM13v2.0.fasta        unzipped fasta
  # CHM13v2.0.fasta.amb    bwa index
  # CHM13v2.0.fasta.ann    bwa index
  # CHM13v2.0.fasta.bwt    bwa index
  # CHM13v2.0.fasta.fai    samtools fai
  # CHM13v2.0.fasta.pac    bwa index
  # CHM13v2.0.fasta.sa     bwa index
  # CHM13v2.0.fasta.sizes  TSV of chromosmes and lengths

  input {
    String name
    File fasta
    RunEnv runenv
  }

  String fasta_bn = "${name}.fasta"
  command <<<
    if [[ "~{fasta}" =~ \.gz$ ]]; then
        gunzip -c ~{fasta} > ~{fasta_bn}
    else
        ln ~{fasta} ~{fasta_bn}
    fi
    samtools dict ~{fasta_bn} -o ~{name}.dict
    samtools faidx ~{fasta_bn} -o ~{fasta_bn}.fai
    cut -f1 ~{fasta_bn}.fai > ~{fasta_bn}.sizes
    bwa index -p ~{fasta_bn} ~{fasta_bn}
    tar cvvf ~{name}.tar ~{name + ".*"}
  >>>

  runtime {
    docker: runenv.docker
    cpu : runenv.cpu
    memory : "~{runenv.memory} GB"
  }

  output {
    File idx = "${name}.tar"
  }
}

task run_untar_idx {
  input {
    File idx
    RunEnv runenv
  }

  command <<<
    mkdir ref
    cd ref
    tar -xvvf ~{idx}
    find . -type f -exec touch {} \;
    find -name \*.sizes | xargs -I% cut -f1 % > chromosomes
  >>>

  output {
    Directory path = "ref"
    File dict = glob("ref/*.dict")[0]
    File fasta = glob("ref/*.fasta")[0]
    File fai = glob("ref/*.fasta.fai")[0]
    File sizes = glob("ref/*.sizes")[0]
    File amb = glob("ref/*.amb")[0]
    File ann = glob("ref/*.ann")[0]
    File bwt = glob("ref/*.bwt")[0]
    File pac = glob("ref/*.pac")[0]
    File sa = glob("ref/*.sa")[0]
    Array[String] chromosomes = read_lines("ref/chromosomes")
  }

  runtime {
    docker: runenv.docker
    cpu : runenv.cpu
    memory : runenv.memory +" GB"
  }
}
