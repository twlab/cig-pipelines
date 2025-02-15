version development

import "../../../structs/runenv.wdl"

task run_verifybamid {
  input {
    File bam
    File ref_fasta
    File ref_fai
    Directory resource
    String svdprefix
    RunEnv runenv
  }

  command <<<
    ln ~{ref_fasta} ./
    ln ~{ref_fai} ./
    verifybamid2 \
      --SVDPrefix ~{resource}/~{svdprefix}
      --Reference ~{ref_fasta} \
      --BamFile ~{bam}
  >>>

  output {
    File ancestry = glob("result.Ancestry")[0]
    File selfsm = glob("result.selfSM")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "${runenv.memory} GB"
    #disks: runenv.disks
  }
}
