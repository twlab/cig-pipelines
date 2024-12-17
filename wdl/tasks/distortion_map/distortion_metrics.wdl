version development

import "../../structs/runenv.wdl"

task calculate_distortion_metrics {
  input {
    File normalized_aligned_reference_matrix
    File normalized_lifted_aligned_source_matrix
    File interval_mapping
    RunEnv runenv
  }

  # --normalized_aligned_reference_matrix       Path to the normalized aligned in reference matrix (.mtx format).
  # --normalized_lifted_aligned_source_matrix   Path to the normalized lifted aligned in source matrix (.mtx format).
  # --intervals_file INTERVALS_FILE             Path to the interval mapping file (TSV format).
  # --output_file OUTPUT_FILE                   Path to the output file (TSV format).
  String output_fn = "distortion_metrics.tsv"
  command <<<
    /apps/scritps/calculate_distortion_metrics.py \
    --normalized_aligned_reference_matrix ~{normalized_aligned_reference_matrix} \
    --normalized_lifted_aligned_source_matrix ~{normalized_lifted_aligned_source_matrix} \
    --intervals_file ~{interval_mapping} \
    --output_file 
  >>>

  output {
    File metrics = glob("~{output_fn}")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
    #disks: runenv.disks
  }
}
