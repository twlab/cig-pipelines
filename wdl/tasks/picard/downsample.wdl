version development

import "../../structs/runenv.wdl"

task run_downsample {
  input {
    File bam
    String probability
    String params = ""
    RunEnv runenv
  }

  String output_bam = basename("~{bam}")
  String metrics_file = "~{output_bam}.metrics"
  command <<<
    java -Xmx~{runenv.memory - 1}g -jar /usr/picard/picard.jar DownsampleSam \
      --INPUT ~{bam} \
      --OUTPUT ~{output_bam} \
      --PROBABILITY ~{probability} \
      --METRICS_FILE ~{metrics_file} \
      ~{params}
  >>>

  output {
    File output_bam = "~{output_bam}"
    File? output_bai = "~{output_bam}.bai"
    File metrics = "~{metrics_file}"
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks: runenv.disks
  }
}

task run_calculate_probability {
  input {
    File samtools_stat_file
    Int genome_size_gb
    Int coverage
    RunEnv runenv
  }

  command <<<
    set -e
    total_bases=$(<~{samtools_stat_file} grep "total length" | awk -F"\t" '{print $3}')
    total_reads=$(<~{samtools_stat_file} grep "total sequences" | awk -F"\t" '{print $3}')
    read_avg_length=$(<~{samtools_stat_file} grep "average length" | awk -F"\t" '{print $3}')

    genome_size=$((~{genome_size_gb} * 1000000000))
    bases_needed=$((~{coverage} * ${genome_size}))
    reads_needed=$(perl -e "print(int(~{coverage} * ${genome_size} / ${read_avg_length}))")
    current_coverage=$(perl -e "printf('%f', (${total_bases} / ${genome_size}))")

    printf "Genome size:          %i\n" "${genome_size}"
    printf "Coverage:             %i\n" "~{coverage}"
    printf "Bases:                %i\n" "${total_bases}"
    printf "Reads:                %i\n" "${total_reads}"
    printf "Average read length:  %i\n" "${read_avg_length}"
    printf "Bases needed:         %i\n" "${bases_needed}"
    printf "Reads needed:         %i\n" "${reads_needed}"
    printf "Current coverage:     %.1fX\n" "${current_coverage}"
    if test "${total_bases}" -lt "${bases_needed}"; then
      echo "Needed bases exceeds total bases. Probability is 1."
      echo 1 > probability
    else
      probability=$(perl -e "printf('%.4f', (${reads_needed} / ${total_reads}))")
      printf "Probability:  ${probability}\n"
      printf "${probability}" > probability
    fi
    printf "Wrote probability to file: probability\n"
    sleep 30s
  >>>

  output {
    Float probability = read_float("probability")
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    #disks : select_first([runenv.disks,"local-disk 100 SSD"])
  }
}
