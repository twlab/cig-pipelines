version development

import "../structs/runenv.wdl"

task load_db {
  input {
    File source_positions
    File lifted_source
    File aligned_ref
    File aligned_source
    File lifted_aligned_source
    RunEnv runenv
  }

  # --db                      Path to the SQLite database file
  # --source_positions        Path to the source positions BED file
  # --lifted_source           Path to the lifted source positions BED file
  # --aligned_ref             Path to the aligned in reference positions BED file
  # --aligned_source          Path to the aligned in source positions BED file
  # --lifted_aligned_source   Path to the lifted aligned in source positions BED file
  command <<<
    Integrate_simulation_data_standardized_into_SQL_db.py \
      --db dm.db \
      --source_positions ~{source_positions} \
      --lifted_source ~{lifted_source} \
      --aligned_ref ~{aligned_ref} \
      --aligned_source ~{aligned_source} \
      --lifted_aligned_source ~{lifted_aligned_source}
  >>>

  output {
    File db = glob("dm.db")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
    #disks: runenv.disks
  }
}
