version development

import "../../structs/runenv.wdl"

task extract_haplotype_reads_from_reference_aligned_cram {
  meta {
    author: "Eddie Belter"
    version: "0.1"
    description: "Extract Haplotype Aligned Reads from a Reference Aligned CRAM."
  }

  input {
    Array[String] haplotype_names
    Array[File] haplotype_reads_fofs
    File crams_and_refs_tsv
    RunEnv runenv
  }

  String hap1_output_cram = "out/~{haplotype_names[0]}.cov.cram"
  String hap2_output_cram = "out/~{haplotype_names[1]}.cov.cram"
  command <<<
    mkdir out/
    printf "%s\t%s\n" "~{haplotype_reads_fofs[0]}" "~{hap1_output_cram}" > haplotypes_and_sources.tsv
    printf "%s\t%s\n" "~{haplotype_reads_fofs[1]}" "~{hap2_output_cram}" >> haplotypes_and_sources.tsv
    cat ~{crams_and_refs_tsv} >> haplotypes_and_sources.tsv
    phase-coverage haplotypes_and_sources.tsv -d out/
  >>>

  output {
    File hap1_cram = hap1_output_cram
    File hap2_cram = hap2_output_cram
    File stats = "out/phase-coverage.stats.yaml"
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
  }
}
