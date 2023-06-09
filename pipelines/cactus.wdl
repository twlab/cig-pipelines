version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/cactus.wdl"

workflow pangenome_graph {
    input {
        File sequences         # tsv of sequences [fasta] with special header
        String reference_name  # name of the sequence to be used as the reference
        String out_name     # base name for outputs
        String docker = "ebelter/cactus:2.5.0-20.04"
        Int cpu = 4
        Int memory = 16
        Int disks = 100
    }

    RunEnv runenv = {
      "docker": docker,
      "cpu": cpu,
      "memory": memory,
      "disks": disks,
    }

    call cactus.minigraph { input:
        sequences=sequences,
        reference_name=reference_name,
        out_name=out_name,
        runenv=runenv,
    }

    call cactus.graphmap { input:
        sequences=sequences,
        gfa=minigraph.gfa,
        reference_name=reference_name,
        out_name=out_name,
        runenv=runenv,
    }

    call cactus.graphmapsplit { input:
        sequences=graphmap.updated_sequences,
        fasta=graphmap.fasta,
        gfa=minigraph.gfa,
        paf=graphmap.paf,
        reference_name=reference_name,
        out_name=out_name,
        runenv=runenv,
    }

    call cactus.align { input:
        chroms=graphmapsplit.chroms,
        fasta=graphmap.fasta,
        reference_name=reference_name,
        out_name=out_name,
        runenv=runenv,
    }

    call cactus.minigraphjoin { input:
        alignments=align.alignments,
        reference_name=reference_name,
        out_name=out_name,
        runenv=runenv,
    }

    output {
        File dist = minigraphjoin.dist
        File gbz = minigraphjoin.gbz
        File gfa = minigraphjoin.gfa
        File hal = minigraphjoin.hal
        File min = minigraphjoin.min
        File stats = minigraphjoin.stats
        File raw_vcf = minigraphjoin.raw_vcf
        File raw_vcf_tbi = minigraphjoin.raw_vcf_tbi
        File vcf = minigraphjoin.vcf
        File vcf_tbi = minigraphjoin.vcf_tbi
    }
}
