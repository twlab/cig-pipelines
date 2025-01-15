version development

import "../../structs/runenv.wdl"

task generate {
  input {
    File db
    File reference_intervals
    RunEnv runenv
  }

  # --db                                       Path to the SQLite database.
  # --threads                                  Number of threads to use.
  # --reference_intervals_file                 Path to the file containing reference intervals.
  # --aligned_reference_count_matrix_file      Path to save the aligned in reference count matrix in .mtx format.
  # --lifted_aligned_source_count_matrix_file  Path to save the lifted aligned in source count matrix in .mtx format.
  # --interval_mapping_file                    Path to save the interval to index mapping file.
  command <<<
    /apps/scripts/generate_count_matrices.py \
      --db ~{db} \
      --threads ~{runenv.cpu} \
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

task normalize {
  input {
    File matrix
    RunEnv runenv
  }

  # --input   Path to the input merged count matrix
  # --output  Path to save the normalized count matrix
  String output_file = basename(matrix, ".mtx") + ".normalized.mtx"
  command <<<
    /apps/scripts/normalize_count_matrices.py \
      --input ~{matrix} \
      --output ~{output_file}
  >>>

  output {
    File normalized_matrix = glob(output_file)[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
    #disks: runenv.disks
  }
}

task merge {
  input {
    Array[File] matrices
    String output_fn
    RunEnv runenv
  }

  #--matrix_file_list  Path to the file containing a list of matrix files. (FOF)
  #--output_file       Path to save the merged matrix in .mtx format.
  command <<<
    /apps/scripts/merge_count_matrices.py \
    --matrix_file_list ~{write_lines(matrices)} \
    --output_file ~{output_fn}
  >>>

  output {
    File matrix = glob("~{output_fn}")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
    #disks: runenv.disks
  }
}
