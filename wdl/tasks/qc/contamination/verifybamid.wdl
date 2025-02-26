version development

import "../../../structs/runenv.wdl"

task run_verifybamid {
  input {
    String sample
    File bam
    File bai
    File ref_fasta
    File ref_fai
    File resource
    RunEnv runenv
  }

  String svdprefix = basename(resource, ".tar")
  command <<<
    ln ~{bam} ./
    ln ~{bai} ./
    ln ~{ref_fasta} ./
    ln ~{ref_fai} ./
    mkdir ./resource
    tar xvvf ~{resource} -C resource
    VerifyBamID \
      --SVDPrefix resource/~{svdprefix} \
      --Reference ~{basename(ref_fasta)} \
      --BamFile ~{basename(bam)} | tee ~{sample}.verifybamid2.txt
  >>>

  output {
    File contamination_report = glob("~{sample}.verifybamid2.txt")[0]
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
