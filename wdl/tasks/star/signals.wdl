version development

import "../../structs/runenv.wdl"

# ebelter 2026-06-18
# This task is pulled from https://github.com/ENCODE-DCC/rna-seq-pipeline/blob/dev/rna-seq-pipeline.wdl and calls bam_tosignals in the ENCODE docker. This repo has not been updated for 4+ years.
# Lastest docker info:
# DOCKER: encodedcc/rna-seq-pipeline
# DIGEST: sha256:02569a0dc8ff397af42ebd4f99a2a9df121cc97010f120b51dcce510e7423abd
# PULL:   encodedcc/rna-seq-pipeline@sha256:02569a0dc8ff397af42ebd4f99a2a9df121cc97010f120b51dcce510e7423abd
task run_star_signals {
  input {
    File bam
    String strandedness      # Stranded, Unstranded
    String ref_prefix = "-"  # chr
    RunEnv runenv
  }

  String bam_bn = basename(bam, ".bam")
  command <<<
    STAR --runMode inputAlignmentsFromBAM \
      --inputBAMfile ~{bam} \
      --outWigType bedGraph \
      --outWigStrand ~{strandedness} \
      --outWigReferencesPrefix ~{ref_prefix}
    for f in $(ls *.bg); do
        mv "${f}" "~{bam_bn}.${f}"
    done
  >>>

  output {
    Array[File] bigwigs = glob("*.bg")
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
  }
}
