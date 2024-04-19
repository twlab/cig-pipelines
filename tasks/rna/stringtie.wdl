version development

import "../../structs/runenv.wdl"

# De Novo
# bulk
# stringtie -o ~{denovo_transcipts_fn} -p ~{cpu} -m 100 ~{bam}
# long
# stringtie --rf -L -v -p 2 -c 1.5 -f 0.05 -o ~{denovo_transcipts_fn} ~{bam}
# 
# Quantification
# bulk
# stringtie -o ~{abundance_estimate_fn} -p ~{cpu} -m 100 -e -b ~{sample} -G ~{annotation_gunzip} ~{bam}
# long
# stringtie --rf -L -v -p 2 -c 1.5 -f 0.05 -e -b ~{ballgown_dn} -G ~{annotation_gunzip} -o ~{abundance_estimate_fn} ~{bam}

# Summary
# BULK: -m 100 
# LONG: -L -c 1.5 -f 0.05
# QUANT: -e -G ~{anno}

task run_stringtie_denovo {
  input {
    String sample
    File bam
    String strandedness # forward => -fr reverse => -rf unstranded => ""
    String params
    RunEnv runenv
  }

  String fr_strandedness_param = if strandedness == "forward" then "--rf" else "--fr" # set reverse as param, check unstranded next line
  String strandedness_param = if strandedness == "unstranded" then "" else fr_strandedness_param
  String denovo_transcipts_fn = "~{sample}.stringtie2.denovo.gtf"
  command <<<
    stringtie \
      ~{bam} \
      ~{strandedness_param} \
      -p ~{runenv.cpu} \
      -o ~{denovo_transcipts_fn} \
      ~{params}
  >>>

  output {
    File denovo_transcipts = glob("~{denovo_transcipts_fn}")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks : runenv.disks
  }
}

task run_stringtie_quantification {
  input {
    String sample
    File bam
    String strandedness # forward => -fr reverse => -rf unstranded => ""
    File annotation # unzipped
    String params
    RunEnv runenv
  }

  String fr_strandedness_param = if strandedness == "forward" then "--rf" else "--fr" # set reverse as param, check unstranded next line
  String strandedness_param = if strandedness == "unstranded" then "" else fr_strandedness_param
  String abundance_estimate_fn = "~{sample}.stringtie2.abundance_estimate.gtf"
  command <<<
    stringtie \
      ~{bam} \
      ~{strandedness_param} \
      -e \
      -G ~{annotation} \
      -b ballgown_tables \
      -p ~{runenv.cpu} \
      -o ~{abundance_estimate_fn} \
      ~{params}
  >>>

  output {
    File abundance_estimate = glob("~{abundance_estimate_fn}")[0]
    Directory ballgown_tables = "ballgown_tables"
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks : runenv.disks
  }
}
