version development

import "../../structs/runenv.wdl"

task run_minibwa {
  input {
    String sample
    String library
    String rg_id
    String platform = "ILLUMINA"
    String platform_unit

    Array[File] fastqs
    Array[File] idx_files
    String minibwa_params = ""

    Int sort_cpu = 2
    Int sort_memory = 4
    Int sort_compression_level = 1

    RunEnv runenv
  }

  String bam = "~{sample}.sorted.bam"
  String reference_fasta_bn = basename(idx_files[0])

  # Minibwa uses -t worker threads plus up to two additional I/O threads.
  # Reserve sort_cpu threads plus two potential minibwa I/O threads.
  Int minibwa_cpu = if runenv.cpu > sort_cpu + 2
    then runenv.cpu - sort_cpu - 2
    else 1

  command <<<
    set -euo pipefail

    if [ "~{length(fastqs)}" -ne 2 ]; then
      echo "ERROR: minibwa paired-end mapping requires exactly two FASTQ files." >&2
      exit 1
    fi

    if [ "~{length(idx_files)}" -ne 3 ]; then
      echo "ERROR: Expected FASTA, L2B, and MBW inputs." >&2
      exit 1
    fi

    if [ "~{sort_cpu}" -lt 1 ]; then
      echo "ERROR: sort_cpu must be at least 1." >&2
      exit 1
    fi

    if [ "~{runenv.cpu}" -lt $((~{sort_cpu} + 2)) ]; then
      echo "ERROR: runenv.cpu must provide sort threads plus two minibwa I/O threads." >&2
      exit 1
    fi

    mkdir -p ref

    ln -sf "~{idx_files[0]}" "ref/~{basename(idx_files[0])}"
    ln -sf "~{idx_files[1]}" "ref/~{basename(idx_files[1])}"
    ln -sf "~{idx_files[2]}" "ref/~{basename(idx_files[2])}"

    reference_fasta="ref/~{reference_fasta_bn}"

    for required_file in \
      "${reference_fasta}" \
      "${reference_fasta}.l2b" \
      "${reference_fasta}.mbw" \
      "~{fastqs[0]}" \
      "~{fastqs[1]}"
    do
      if [ ! -s "${required_file}" ]; then
        echo "ERROR: Missing or empty input: ${required_file}" >&2
        exit 1
      fi
    done

    minibwa map \
      -x sr \
      -t ~{minibwa_cpu} \
      -K 320000000 \
      -R '@RG\tID:~{rg_id}\tLB:~{library}\tSM:~{sample}\tPL:~{platform}\tPU:~{platform_unit}' \
      ~{minibwa_params} \
      "${reference_fasta}" \
      "~{fastqs[0]}" \
      "~{fastqs[1]}" \
    | samtools sort \
        -@ ~{sort_cpu} \
        -m "~{sort_memory}"G \
        -l ~{sort_compression_level} \
        -T "~{sample}.sorttmp" \
        -O BAM \
        -o "~{bam}" \
        -

    samtools quickcheck -v "~{bam}"
  >>>

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
  }

  output {
    File bam = "~{bam}"
  }
}
