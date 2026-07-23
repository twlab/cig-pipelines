version development

import "../../structs/runenv.wdl"

task run_build_idx {
  meta {
    description: "Build a minibwa reference index and package it with reference support files"
  }

  input {
    String name
    File fasta
    RunEnv runenv
  }

  String fasta_bn = "~{name}.fasta"

  command <<<
    set -xeuo pipefail

    if [ "~{runenv.cpu}" -lt 4 ]; then
      echo "ERROR: runenv.cpu must provide at least 4 threads for minibwa indexing." >&2
      exit 2
    fi

    mkdir -p ref
    cd ref

    if [[ "~{fasta}" =~ \.gz$ ]]; then
      gunzip -c "~{fasta}" > "~{fasta_bn}"
    else
      # Store a real FASTA in the archive, not a symlink into a Cromwell call dir.
      cp --reflink=auto "~{fasta}" "~{fasta_bn}"
    fi

    if [ ! -s "~{fasta_bn}" ]; then
      echo "ERROR: Reference FASTA is missing or empty: ~{fasta_bn}" >&2
      exit 1
    fi

    samtools dict \
      "~{fasta_bn}" \
      -o "~{name}.dict"

    samtools faidx "~{fasta_bn}"

    awk 'BEGIN {OFS="\t"} {print $1, $2}' \
      "~{fasta_bn}.fai" \
      > "~{fasta_bn}.sizes"

    minibwa index \
      -t ~{runenv.cpu} \
      "~{fasta_bn}"

    for required_file in \
      "~{name}.dict" \
      "~{fasta_bn}" \
      "~{fasta_bn}.fai" \
      "~{fasta_bn}.sizes" \
      "~{fasta_bn}.l2b" \
      "~{fasta_bn}.mbw"
    do
      if [ ! -s "${required_file}" ]; then
        echo "ERROR: Missing or empty required file: ${required_file}" >&2
        exit 1
      fi
    done

    tar \
      vvcf "../~{name}.tar" \
      "~{name}.dict" \
      "~{fasta_bn}" \
      "~{fasta_bn}.fai" \
      "~{fasta_bn}.sizes" \
      "~{fasta_bn}.l2b" \
      "~{fasta_bn}.mbw"

    tar tvvf "../~{name}.tar"
  >>>

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
  }

  output {
    File idx = "~{name}.tar"
    File dict = "ref/~{name}.dict"
    File fai = "ref/~{fasta_bn}.fai"
    File FASTA = "ref/~{fasta_bn}"
    File sizes = "ref/~{fasta_bn}.sizes"
    File l2b = "ref/~{fasta_bn}.l2b"
    File mbw = "ref/~{fasta_bn}.mbw"
  }
}

task run_untar_idx {
  meta {
    description: "Extract and validate a minibwa reference-index archive"
  }

  input {
    File idx
    RunEnv runenv
  }

  command <<<
    set -euo pipefail

    mkdir -p ref
    cd ref

    tar -xf "~{idx}"
    find . -type f -exec touch {} +

    fasta_count=$(find . -maxdepth 1 -type f -name "*.fasta" | wc -l)
    dict_count=$(find . -maxdepth 1 -type f -name "*.dict" | wc -l)
    sizes_count=$(find . -maxdepth 1 -type f -name "*.sizes" | wc -l)
    l2b_count=$(find . -maxdepth 1 -type f -name "*.l2b" | wc -l)
    mbw_count=$(find . -maxdepth 1 -type f -name "*.mbw" | wc -l)

    if [ "${fasta_count}" -ne 1 ] || \
       [ "${dict_count}" -ne 1 ] || \
       [ "${sizes_count}" -ne 1 ] || \
       [ "${l2b_count}" -ne 1 ] || \
       [ "${mbw_count}" -ne 1 ]; then
      echo "ERROR: Reference archive does not contain exactly one FASTA, DICT, sizes, L2B, and MBW file." >&2
      find . -maxdepth 1 -type f -print | sed 's#^\./##' | sort >&2
      exit 1
    fi

    reference_fasta=$(find . -maxdepth 1 -type f -name "*.fasta" -print -quit)
    dictionary=$(find . -maxdepth 1 -type f -name "*.dict" -print -quit)
    sizes_file=$(find . -maxdepth 1 -type f -name "*.sizes" -print -quit)

    for required_file in \
      "${reference_fasta}" \
      "${reference_fasta}.fai" \
      "${reference_fasta}.sizes" \
      "${reference_fasta}.l2b" \
      "${reference_fasta}.mbw" \
      "${dictionary}"
    do
      if [ ! -s "${required_file}" ]; then
        echo "ERROR: Missing or empty reference component: ${required_file}" >&2
        exit 1
      fi
    done

    cut -f1 "${sizes_file}" > chromosomes
    test -s chromosomes

    find . -maxdepth 1 -type f -print | sed 's#^\./##' | sort
  >>>

  output {
    # Retained for compatibility with directory-oriented Cromwell workflows.
    Directory path = "ref"

    File dict = glob("ref/*.dict")[0]
    File fasta = glob("ref/*.fasta")[0]
    File fai = glob("ref/*.fasta.fai")[0]
    File sizes = glob("ref/*.sizes")[0]
    File l2b = glob("ref/*.l2b")[0]
    File mbw = glob("ref/*.mbw")[0]
    Array[String] chromosomes = read_lines("ref/chromosomes")
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
  }
}
