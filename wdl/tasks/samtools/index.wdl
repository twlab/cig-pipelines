version development

import "../../structs/runenv.wdl"

task run_index {
  input {
    File bam
    RunEnv runenv
  }

  String bai = "~{basename(bam)}.bai"
  command <<<
    ln ~{bam} ~{basename(bam)}
    samtools index -b -@~{runenv.cpu} -o ~{bai} ~{basename(bam)}
  >>>

  output {
    File bai = glob(bai)[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
  }
}

task new_run_index {
  input {
    File sam_file
    File? reference
    RunEnv runenv
  }

  String sam_bn = basename(sam_file)
  String ext = sub(sam_bn, "^.*\\.", "")
  String idx_ext = sub(ext, "am$", "ai")
  String idx = "~{sam_bn}.~{idx_ext}"

  command <<<
    set -ex
    #sam_bn=$(basename "~{sam_file}")
    #ln "~{sam_file}" "${sam_bn}"
    #ext="${sam_bn##*.}"
    #idx="${sam_bn}.${ext/%am/ai}"
    #samtools index -b -@ ~{runenv.cpu} -o "${idx}" ~{if (defined(reference)) then "--reference ~{reference}" else ""} ~{basename(sam_file)}
    ln "~{sam_file}" "~{sam_bn}"
    samtools index -b -@ ~{runenv.cpu} -o "~{idx}" ~{if (defined(reference)) then "--reference ~{reference}" else ""} ~{sam_bn}
  >>>

  output {
    File idx = idx
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
  }
}
