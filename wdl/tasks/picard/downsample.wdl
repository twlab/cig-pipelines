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
    File? output_bai = glob("*.bai")[0]
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
    set -ex
    genome_size=$((~{genome_size_gb} * 1000000000))
    bases_needed=$((~{coverage} * ${genome_size}))
    total_bases=$(<~{samtools_stat_file} grep "total length" | awk -F"\t" '{print $3}')
    total_reads=$(<~{samtools_stat_file} grep "total sequences" | awk -F"\t" '{print $3}')
    read_avg_length=$(<~{samtools_stat_file} grep "average length" | awk -F"\t" '{print $3}')
    reads_needed=$(perl -e "print(int(~{coverage} * ${genome_size} / ${read_avg_length}))")

    printf "Genome size:   %i\n" "${genome_size}"
    printf "Coverage:      %i\n" "~{coverage}"
    printf "Total bases:   %i\n" "${total_bases}"
    printf "Bases needed:  %i\n" "${bases_needed}"
    printf "Total reads:   %i\n" "${total_reads}"
    printf "Reads needed:  %i\n" "${reads_needed}"
    if test "${total_bases}" -lt "${bases_needed}"; then
      echo "Needed bases exceeds total bases. Probability is 1."
      echo 1 > probability
    else
      probability=$(perl -e "print(int(${reads_needed} / ${total_reads}))")
      printf "Probability:  ${probability}\n"
      printf "${probability}" > probability
    fi
    printf "Wrote probability to file: probability\n"
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
