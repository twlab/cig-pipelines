version development

import "../../structs/runenv.wdl"

task run_cpg {
  input {
    File alignments       # BAM or CRAM
    File alignments_idx   # BAI or CRAI
    File reference_fasta  # for CRAM
    String params = "" 
    RunEnv runenv
  }

  # v3.0.0 - no model needed
  # --mode [default: model] [possible values: count, model]
  # --min-coverage[default: 4]
  # --min-mapq [default: 1]
  # 
  # DEPRECATED v2.3.2 - model location required
  # --model /opt/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite

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
    File bed_tbi = "~{output_prefix}.combined.bed.gz.tbi"
    File bigwig = "~{output_prefix}.combined.bw"
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks : runenv.disks
  }
}
