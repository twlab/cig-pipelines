version development

import "../structs/runenv.wdl"

task run_samblaster {
  input {
    File bam
    RunEnv runenv
  }

  String output_bam_bn = basename(bam, ".bam")
  String output_bam = "~{output_bam_bn}.dedup.bam"
  # -i --input           FILE Input sam file [stdin].
  # -o --output          FILE Output sam file for all input alignments [stdout].
  # -M                        Run in compatibility mode; both 0x100 and 0x800 are considered chimeric. Similar to BWA MEM -M option.
  # -a --acceptDupMarks       Accept duplicate marks already in input file instead of looking for duplicates in the input.
  # -e --excludeDups          Exclude reads marked as duplicates from discordant, splitter, and/or unmapped file.
  # -r --removeDups           Remove duplicates reads from all output files. (Implies --excludeDups).
  # --addMateTags          Add MC and MQ tags to all output paired-end SAM lines.
  # --ignoreUnmated        Suppress abort on unmated alignments. Use only when sure input is read-id grouped,
  # -q --quiet                Output fewer statistics.
  command <<<
    samblaster -M --input ~{bam} --output ~{output_bam}
  >>>

  output {
    File dedup_bam = glob("~{output_bam}")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
  }
}
