version development

import "../../structs/runenv.wdl"

task run_align {
  input {
    File bam
    File reference_mmi
    RunEnv runenv
  }

  String output_bam = basename(bam)
  Int threads = runenv.cpu - 4
  Int sort_threads = runenv.cpu - threads
  command <<<
    set -eo pipefail
    mkdir tmpsort/
    pbmm2 align --preset CCS -j ~{threads} --log-level INFO ~{reference_mmi} ~{bam} | samtools sort -@ ~{sort_threads} -T tmpsort/sorted.nnnn.bam -O bam -o ~{output_bam}##idx##~{output_bam}.bai --write-index
    rm -rf tmpsort/
  >>>

  output {
    File aligned_bam = glob("~{output_bam}")[0]
    File aligned_bai = glob("~{output_bam}.bai")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks : runenv.disks
  }
}
