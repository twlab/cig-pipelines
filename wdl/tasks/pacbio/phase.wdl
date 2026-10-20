version development

import "../../structs/runenv.wdl"

task run_phase_crams {
  meta {
    author: "Eddie Belter"
    version: "1.0"
  }

  input {
    File haplotype1_cram
    File haplotype2_cram
    String haplotype1_reference
    String haplotype2_reference
    RunEnv runenv
  }

  String haplotype1_cram_bn = basename(haplotype1_cram)
  String haplotype2_cram_bn = basename(haplotype2_cram)
  command <<<
    set -e
    touch crams_and_refs.tsv
    ln "~{haplotype1_cram}" "~{haplotype1_cram_bn}"
    ln "~{haplotype2_cram}" "~{haplotype2_cram_bn}"
    printf "%s\t%s\n" "~{haplotype1_cram_bn}" "~{haplotype1_reference}" >> crams_and_refs.tsv
    printf "%s\t%s\n" "~{haplotype2_cram_bn}" "~{haplotype2_reference}" >> crams_and_refs.tsv
    /usr/local/bin/phase-reads crams_and_refs.tsv
  >>>

  output {
    Array[File] phased_crams = glob("*.phased.cram")
    Array[File] yamls = glob("*.yaml")
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
  }
}
