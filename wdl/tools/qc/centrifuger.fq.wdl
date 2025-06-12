version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/qc/contamination/centrifuger.wdl"

workflow centrifuger {
  input {
    Array[File] fastqs
    File centrifuger_db
    File centrifuger_taxon_db
    String centrifuger_cutoff_percentage
    String docker
    Int cpu
    Int memory
  }

  RunEnv runenv_centrifuger = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
    "disks": 20,
  }

  call centrifuger.run_centrifuger { input:
    fastqs=fastqs,
    centrifuger_db=centrifuger_db,
    centrifuger_taxon_db=centrifuger_taxon_db,
    cutoff_percentage=centrifuger_cutoff_percentage,
    runenv=runenv_centrifuger,
  }
}
