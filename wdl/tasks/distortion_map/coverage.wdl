version development

import "../../structs/runenv.wdl"

task generate_simulated_coverage {
  input {
    File db
    File simulated_intervals
    Int batch_size
    RunEnv runenv
  }

  # --db DB           Path to the SQLite database.
  # --intervals_file  Path to the intervals file.
  #--output_file      Output file for the simulated coverage in BED format.
  # --batch_size      Batch size for processing source positions.
  command <<<
    /scratch/hllab/Juan/JuanMacias_General_Code/DistortionMapping/calculate_simulated_source_coverage.py \
      --db ~{db} \
      --intervals_file ~{simulated_intervals} \
      --batch_size ~{batch_size} \
      --output_file simulated_coverage.bed
  >>>

  output {
    File coverage = glob("simulated_coverage.bed")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
    #disks: runenv.disks
  }
}

task generate_simulated_no_lift_over_coverage {
  input {
    File db
    File simulated_intervals
    RunEnv runenv
  }

  # --db             Path to the SQLite database.
  # --intervals_file Path to the intervals file.
  # --output_file    Output file for the simulated no-lift-over coverage in BED format.
  command <<<
    /apps/scripts/calculate_no_lift_over_coverage.py \
      --db ~{db} \
      --intervals_file ~{simulated_intervals} \
      --output_file simulated_no_lift_over_coverage.bed
  >>>

  output {
    File coverage = glob("simulated_no_lift_over_coverage.bed")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
    #disks: runenv.disks
  }
}
