version development

import "../../structs/runenv.wdl"

task generate_count_matrices {
  input {
    File db
    File reference_intervals
    RunEnv runenv
  }

  # --db                                       Path to the SQLite database.
  # --reference_intervals_file                 Path to the file containing reference intervals.
  # --aligned_reference_count_matrix_file      Path to save the aligned in reference count matrix in .mtx format.
  # --lifted_aligned_source_count_matrix_file  Path to save the lifted aligned in source count matrix in .mtx format.
  # --interval_mapping_file                    Path to save the interval to index mapping file.
  command <<<
    /apps/scripts/generate_count_matrices.py \
      --db ~{db} \
      --reference_intervals_file ~{reference_intervals} \
      --aligned_reference_count_matrix_file aligned_reference_count.mtx \
      --lifted_aligned_source_count_matrix_file lifted_aligned_source_count.mtx \
      --interval_mapping_file interval_mapping.tsv
  >>>

  output {
    File aligned_reference_count_matrix = glob("aligned_reference_count.mtx")[0]
    File lifted_aligned_source_count_matrix = glob("lifted_aligned_source_count.mtx")[0]
    File interval_mapping = glob("interval_mapping.tsv")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
    #disks: runenv.disks
  }
}
