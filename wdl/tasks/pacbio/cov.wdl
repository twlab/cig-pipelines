version development

import "../../structs/runenv.wdl"

task extract_haplotype_reads_from_reference_aligned_cram {
  meta {
    author: "Eddie Belter"
    version: "0.1"
    description: "Extract Haplotype Aligned Reads from a Reference Aligned CRAM."
  }

  input {
    File haplotype1_aln_cram
    File haplotype1_aln_ref
    File haplotype2_aln_cram
    File haplotype2_aln_ref
    File reference_aln_cram
    File reference_aln_ref
    RunEnv runenv
  }

  command <<<
    mkdir out/
    printf "%s\t%s\n" ~{haplotype1_aln_cram} ~{haplotype1_aln_ref} > crams_and_ref.tsv
    printf "%s\t%s\n" ~{haplotype2_aln_cram} ~{haplotype2_aln_ref} >> crams_and_ref.tsv
    printf "%s\t%s\n" ~{reference_aln_cram} ~{reference_aln_ref} >> crams_and_ref.tsv
    phase-coverage TSV crams_and_refs.tsv -d out/
  >>>

  output {
    File haplotype1_coverage_cram = basename(haplotype1_aln_cram, ".cram") + ".cov.cram"
    File haplotype2_coverage_cram = basename(haplotype2_aln_cram, ".cram") + ".cov.cram"
    File stats = "out/phase-coverage.stats.yaml"
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
  }
}
