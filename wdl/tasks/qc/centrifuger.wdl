version development

import "../../structs/runenv.wdl"

task run_centrifuger {
  input {
    Array[File] fastqs
    Directory centrifuger_db_path
    Directory centrifuger_taxon_db_path
    String cutoff_percentage
    RunEnv runenv
  }

  command <<<
    set -x
    centrifuger_db_bn=$(find ~{centrifuger_db_path} -name '*.cfr' | xargs -I% basename % | awk -F. '{print $1}' | uniq)
    centrifuger_db="~{centrifuger_db_path}/${centrifuger_db_bn}"
    centrifuger -x "${centrifuger_db}" -1 ~{fastqs[0]} -2 ~{fastqs[1]} -t ~{runenv.cpu} > centrifuger.output.tsv
    summarize_centrifuger.py --centrifuger_output centrifuger.output.tsv --taxon_db ~{centrifuger_taxon_db_path} --summary_output centrifuger.summary --threads ~{runenv.cpu} --cutoff_percentage ~{cutoff_percentage} --csv
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
