version development

import "../../structs/runenv.wdl"

task run_bam_to_cram {
  input {
    File bam
    File ref_fasta
    RunEnv runenv
  }

  String output_fn = basename(bam, ".bam") + ".cram"
  command <<<
    samtools view -@ ~{runenv.cpu} -h -o ~{output_fn} -C -T ~{ref_fasta} ~{bam}
  >>>

  output {
    File cram = "~{output_fn}"
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
  }
}
