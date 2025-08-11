version development

import "../../structs/runenv.wdl"

task run_align {
  input {
    String sample
    File bam
    File reference_mmi
    String params
    RunEnv runenv
  }

  String output_bam = "~{sample}.bam"
  Int threads = runenv.cpu - 4
  Int sort_threads = runenv.cpu - threads
  command <<<
    set -eo pipefail
    mkdir tmpsort/
    pbmm2 align ~{params} --sample ~{sample} -j ~{threads} --log-level INFO ~{reference_mmi} ~{bam} | samtools sort -@ ~{sort_threads} -T tmpsort/sorted.nnnn.bam -O bam -o ~{output_bam}##idx##~{output_bam}.bai --write-index
    rm -rf tmpsort/
  >>>

  output {
    File aligned_bam = "~{output_bam}"
    File aligned_bai = "~{output_bam}.bai"
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks : runenv.disks
  }
}

task run_align_output_cram {
  input {
    String sample
    File bam
    File reference_mmi
    File reference_fasta
    String params
    RunEnv runenv
  }

  String output_cram = "~{sample}.cram"
  Int threads = runenv.cpu - 4
  Int sort_threads = runenv.cpu - threads
  command <<<
    set -eo pipefail
    mkdir tmpsort/
    pbmm2 align ~{params} --sample ~{sample} -j ~{threads} --log-level INFO ~{reference_mmi} ~{bam} | samtools sort -@ ~{sort_threads} -T tmpsort/sorted.nnnn.bam -O CRAM --reference ~{reference_fasta} -o ~{output_cram}
    rm -rf tmpsort/
  >>>

  output {
    File aligned_cram = "~{output_cram}"
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks : runenv.disks
  }
}
