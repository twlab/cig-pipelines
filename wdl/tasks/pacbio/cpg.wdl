version development

import "../../structs/runenv.wdl"

task run_cpg_bam {
  input {
    File bam
    File bai
    String model = "/opt/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite"
    String params = "" # --min-coverage[default: 4] --min-mapq [default: 1]
    RunEnv runenv
  }

  String output_prefix = basename(bam, ".bam")
  command <<<
    ln -s ~{bam} .
    ln -s ~{bai} .
    aligned_bam_to_cpg_scores \
      --model ~{model} \
      --bam ~{basename(bam)} \
      --output-prefix ~{output_prefix} \
      --threads ~{runenv.cpu} \
      ~{params}
  >>>

  output {
    File bed = glob("~{output_prefix}*.bed")[0]
    File bigwig = glob("~{output_prefix}*.bw")[0]
    # also has log file
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks : runenv.disks
  }
}

task run_cpg {
  input {
    File alignments       # BAM or CRAM
    File alignments_idx   # BAI or CRAI
    File reference_fasta  # for CRAM
    String params = "" 
    RunEnv runenv
  }

  # Params
  # --mode [default: model] [possible values: count, model]
  # --model /opt/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite
  # --min-coverage[default: 4]
  # --min-mapq [default: 1]

  #String output_prefix1 = basename(alignments, ".bam")
  #String output_prefix = basename(alignments, ".cram")
  String output_prefix = basename(basename(alignments, ".bam"), ".cram")
  command <<<
    ln -s ~{alignments} .
    ln -s ~{alignments_idx} .
    aligned_bam_to_cpg_scores \
      --bam ~{basename(alignments)} \
      ~{if defined(reference_fasta) then "--ref ~{reference_fasta}" else ""} \
      --output-prefix ~{output_prefix} \
      --threads ~{runenv.cpu} \
      ~{params}
  >>>

  output {
    File bed = "~{output_prefix}.combined.bed"
    File bigwig = "~{output_prefix}.combined.bw"
    #File bed = glob("~{output_prefix}*.bed")[0]
    #File bigwig = glob("~{output_prefix}*.bw")[0]
    # also has log file
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks : runenv.disks
  }
}
