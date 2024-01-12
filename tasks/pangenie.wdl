version development

import "../structs/runenv.wdl"

# version: v3.0.1
# usage:
# PanGenie [options] -f <index-prefix> -i <reads.fa/fq> -o <outfile-prefix>
# PanGenie [options] -i <reads.fa/fq> -r <reference.fa> -v <variants.vcf> -o <outfile-prefix>
# 
# options:
#         -a VAL  sample subsets of paths of this size (default: 0).
#         -c      count all read kmers instead of only those located in graph
#         -e VAL  size of hash used by jellyfish (default: 3000000000).
#         -f VAL  Filename prefix of files computed by PanGenie-index (i.e. option -o used with PanGenie-index).
#         -g      run genotyping (Forward backward algorithm, default behaviour)
#         -i VAL  sequencing reads in FASTA/FASTQ format or Jellyfish database in jf format. NOTE: INPUT FASTA/Q FILE MUST NOT BE COMPRESSED.
#         -j VAL  number of threads to use for kmer-counting (default: 1).
#         -k VAL  kmer size (default: 31).
#         -o VAL  prefix of the output files. NOTE: the given path must not include non-existent folders (default: result).
#         -p      run phasing (Viterbi algorithm). Experimental feature
#         -r VAL  reference genome in FASTA format. NOTE: INPUT FASTA FILE MUST NOT BE COMPRESSED.
#         -s VAL  name of the sample (will be used in the output VCFs) (default: sample).
#         -t VAL  number of threads to use for core algorithm. Largest number of threads possible is the number of chromosomes given in the VCF (default: 1).
#         -u      output genotype ./. for variants not covered by any unique kmers
#         -v VAL  variants in VCF format. NOTE: INPUT VCF FILE MUST NOT BE COMPRESSED.
# 
# 
task run_genotyper {
  input {
     File fastq
     Directory index
     String sample
     RunEnv runenv
  }

  #String index_name = basename(index)
  #String index_subd = "~{index}/~{index_name}"
  command <<<
    set -x
    index_name=$(basename ~{index})
    index_subd="~{index}/${index_name}"
    PanGenie -i ~{fastq} -f ${index_subd} -s ~{sample} -o ~{sample} -t ~{runenv.cpu} -j ~{runenv.cpu}
  >>>

  output {
    File vcf = glob("~{sample}*.vcf")[0]
    File histo = glob("~{sample}*.histo")[0]
    File path_segments = glob("~{sample}*.fasta")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "${runenv.memory} GB"
    #disks: runenv.disks
  }
}

task run_index {
  # 8 CPU >64G
  input {
     File ref_fasta
     File graph_vcf
     String name
     RunEnv runenv
  }

  String gunzip_vcf = sub(basename(graph_vcf), ".gz", "")
  String gunzip_ref = sub(basename(ref_fasta), ".gz", "")
  command <<<
    set -ex
    gunzip -c ~{graph_vcf} > gunzip_vcf
    gunzip -c ~{ref_fasta} > gunzip_ref
    PanGenie-index -o ~{name} -r ~{gunzip_ref} -v ~{gunzip_vcf} -t ~{runenv.cpu}
  >>>

  output {
    Directory index = "~{name}"
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "${runenv.memory} GB"
    #disks: runenv.disks
  }
}
