version development

import "../../structs/runenv.wdl"

task create_intervals {
  input {
    File db
    File reference_sizes
    Int window_length
    Int window_stride
    RunEnv runenv
  }

  # --db                        Path to the SQLite database file
  # --input_source_sizes_file   Path to a pre-generated source chromosome sizes file
  # --reference_sizes_file      Path to the reference chromosome sizes file
  # --window_length             Length of the sliding window (default: 24000)
  # --window_stride             Stride of the sliding window (default: 24000)
  # --reference_intervals_file  Output file for reference intervals
  # --simulated_intervals_file  Output file for simulated intervals
  # --source_sizes_file         Output file for source chromosome sizes
  command <<<
    /apps/scripts/prepare_interval_files.py \
      --db ~{db} \
      --reference_sizes_file ~{reference_sizes} \
      --window_length ~{window_length} \
      --window_stride ~{window_stride} \
      --reference_intervals_file reference_intervals.tsv \
      --simulated_intervals_file simulated_intervals.tsv \
      --source_sizes_file source_sizes.tsv
  >>>

  output {
    File reference_intervals = glob("reference_intervals.tsv")[0]
    File simulated_intervals = glob("simulated_intervals.tsv")[0]
    File source_sizes = glob("source_sizes.tsv")[0]
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
    Array[File] intervals
    RunEnv runenv
  }

  command <<<
    for f in ~{sep=' ' intervals}; do
      echo "${f}" >> intervals.fof
    done
    /apps/scripts/merge_interval_files.py \
    intervals.fof \
    merged_intervals.tsv
  >>>

  output {
    File merged_intervals = glob("merged_intervals.tsv")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
    #disks: runenv.disks
  }
}
