version development

import "../../structs/runenv.wdl"

task run_genomecov {
  meta {
    author: "Eddie Belter"
    version: "0.1"
    description: "Run Bedtools Genome Coverage on a SAM/BAM/CRAM"
  }

  input {
    File sam
    File ref
    File output_bg
    RunEnv runenv
  }

  command <<<
    ~{if (defined(ref)) then "CRAM_REFERENCE=~{ref} " else ""} bedtools genomecov -ibam ~{sam} -bga > ~{output_bg}
  >>>

  output {
    File bg = output_bg
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
  }
}

task run_extract_and_genomecov {
  meta {
    author: "Eddie Belter"
    version: "0.1"
    description: "Extract Reads From a SAM/BAM/CRAM then Run Bedtools Genome Coverage"
  }

  input {
    File sam
    File? ref
    File reads_fof
    RunEnv runenv
  }

  String output_file = basename(sam, ".sam") + ".bg"
  command <<<
    samtools view -h --qname-file ~{reads_fof} ~{if (defined(ref)) then "--reference ~{ref}" else ""} ~{sam} | samtools sort -@ ~{runenv.cpu - 1} -O bam - | bedtools genomecov -ibam - -bga > ~{output_file}
  >>>

  output {
    File bg = output_file
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
  }
}
