version development

import "../../structs/runenv.wdl"

task run_slop {
  input {
    File bed_file
    File reference_fai
    Int bases
    RunEnv runenv
  }

  # -b	Increase the BED/GFF/VCF entry by the same number base pairs in each direction. Integer.
  # -l	The number of base pairs to subtract from the start coordinate. Integer.
  # -r	The number of base pairs to add to the end coordinate. Integer.
  # -s	Define -l and -r based on strand. For example. if used, -l 500 for a negative-stranded feature, it will add 500 bp to the end coordinate.
  # -pct	Define -l and -r as a fraction of the feature’s length. E.g. if used on a 1000bp feature, -l 0.50, will add 500 bp “upstream”. Default = false.
  # -header	Print the header from the input file prior to results.
  String output_fn = sub(basename("~{bed_file}"), ".bed$", "") + ".slopped.bed"
  command <<<
    set -eu
    set -o xtrace
    bedtools slop -i "~{bed_file}" -g "~{reference_fai}" -b "~{bases}" > "~{output_fn}"
  >>>

  output {
      File slopped_bed_file = glob("~{output_fn}")[0]
  }

  runtime {
      docker: runenv.docker
      cpu: runenv.cpu
      memory: runenv.memory + " GB"
      #disks: "local-disk " + in_call_disk + " SSD"
  }
}
