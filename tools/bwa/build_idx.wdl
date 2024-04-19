version development

# Create an BWA index TAR with support files.
# CHM13v2.0.dict        samtools dict
# CHM13v2.0.fasta       unzipped fasta
# CHM13v2.0.fasta.amb   bwa index
# CHM13v2.0.fasta.ann   bwa index
# CHM13v2.0.fasta.bwt   bwa index
# CHM13v2.0.fasta.fai   samtools fai
# CHM13v2.0.fasta.pac   bwa index
# CHM13v2.0.fasta.sa    bwa index

import "wdl/structs/runenv.wdl"
import "wdl/tasks/bwa/idx.wdl"

workflow bwa_build_idx {
    meta {
        author: "Eddie Belter"
        version: "0.1"
        description: "Build BWA index from FASTA reference"
    }

    input {
        String name
        File fasta_gz
        String docker = "ebelter/bwa:0.7.17"
        Int cpu = 2
        Int memory = 8
    }

    RunEnv runenv = {
      "docker": docker,
      "cpu": cpu,
      "memory": memory,
      "disks": 20,
    }

    call idx.run_build_idx { input:
        name=name,
        fasta_gz=fasta_gz,
        runenv=runenv
    }

    output {
        File idx = run_build_idx.idx
    }
}
