version development

import "../../../structs/runenv.wdl"

task run_centrifuger {
  input {
    Array[File] fastqs
    File centrifuger_db
    File centrifuger_taxon_db
    String cutoff_percentage
    RunEnv runenv
  }
 
  command <<<
    set -ex
    mkdir db
    tar xvvf ~{centrifuger_db} -C db/
    centrifuger_db_bn=$(find db -name '*.cfr' | xargs -I% basename % | awk -F. '{print $1}' | uniq)
    centrifuger_db_prefix="db/${centrifuger_db_bn}"
    centrifuger -x "${centrifuger_db_prefix}" -1 ~{fastqs[0]} -2 ~{fastqs[1]} -t ~{runenv.cpu} > centrifuger.output.tsv

    mkdir taxon-db
    tar xvvf ~{centrifuger_taxon_db} -C taxon-db/
    summarize_centrifuger.py --centrifuger_output centrifuger.output.tsv --taxon_db taxon-db --summary_output centrifuger.summary --threads ~{runenv.cpu} --cutoff_percentage ~{cutoff_percentage} --csv
  >>>

  output {
    File centrifuger_output_file = glob("centrifuger.output.tsv")[0]
    File summary_output_file = glob("centrifuger.summary.csv")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "${runenv.memory} GB"
    #disks: runenv.disks
  }
}

task run_centrifuger_bam {
  input {
    File bam
    File centrifuger_db
    File centrifuger_taxon_db
    String cutoff_percentage
    RunEnv runenv
  }

  command <<<
    set -x
    mkdir db
    tar xvvf ~{centrifuger_db} -C db/
    centrifuger_db_bn=$(find db -name '*.cfr' | xargs -I% basename % | awk -F. '{print $1}' | uniq)
    centrifuger_db_prefix="db/${centrifuger_db_bn}"
    samtools fastq ~{bam} -@ ~{runenv.cpu} -n -0 R0.fastq -1 R1.fastq -2 R2.fastq
    declare -a fq_params
    if test -s R1.fastq; then
        if test -s R2.fastq; then
            fq_params=("-1" "R1.fastq" "-2" "R2.fastq")
        else
            fq_params=("-u" "R1.fastq")
        fi
    elif test -s R0.fastq; then
        fq_params=("-u" "R0.fastq")
    else
        echo "Could not find FASTQs generated from ~{bam}"
        exit 1
    fi
    centrifuger -x "${centrifuger_db_prefix}" "${fq_params[@]}" -t ~{runenv.cpu} > centrifuger.output.tsv

    mkdir taxon-db
    tar xvvf ~{centrifuger_taxon_db} -C taxon-db/
    summarize_centrifuger.py --centrifuger_output centrifuger.output.tsv --taxon_db taxon-db --summary_output centrifuger.summary --threads ~{runenv.cpu} --cutoff_percentage ~{cutoff_percentage} --csv
  >>>

  output {
    File centrifuger_output_file = glob("centrifuger.output.tsv")[0]
    File summary_output_file = glob("centrifuger.summary.csv")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "${runenv.memory} GB"
    #disks: runenv.disks
  }
}
