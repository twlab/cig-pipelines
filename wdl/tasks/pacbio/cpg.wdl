version development

import "../../structs/runenv.wdl"

task run_cpg {
  input {
    File bam
    File bai
    String model = "/opt/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite"
    RunEnv runenv
  }

  String output_prefix = basename(bam, ".bam") 
  command <<<
    ln -s ~{bam} .
    ln -s ~{bai} .
    aligned_bam_to_cpg_scores \
      --model ~{model} \
      --bam ~{basename(bam)} \
      --output-prefix ~{output_prefix} \
      --threads ~{runenv.cpu}
  >>>

  output {
    File bed = glob("~{output_prefix}*.bed")[0]
    File bigwig = glob("~{output_prefix}*.bw")[0]
    # also has log file
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks : runenv.disks
  }
}
