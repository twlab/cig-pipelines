version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/vcallers/utils.wdl"

workflow validate_vcfs {
  meta {
      author: "Eddie Belter"
      version: "0.1"
      description: "Run a quick validation on a group of VCFs"
  }

  input {
      File vcfs_fof
      String docker = "mgibio/seqtools:latest-noble"
      Int cpu = 2
      Int memory = 10
  }

  RunEnv runenv = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
    "disks": 20,
  }

  Array[Array[String]] samples_and_vcfs = read_tsv(vcfs_fof)
  scatter(sample_and_vcf in samples_and_vcfs) {
    String sample = sample_and_vcf[0]
    File vcf = sample_and_vcf[1]
    call utils.run_quick_validate_vcf as validate { input:
      sample=sample,
      vcf=vcf,
      runenv=runenv
    }
  }

  Array[Array[String]] samples_statuses_and_vcfs = [ validate.sample2, validate.status, validate.vcf2 ]
  call write_report { input:
    samples_statuses_and_vcfs=samples_statuses_and_vcfs,
    runenv=runenv
  }

  output {
  }
}

task write_report {
  input {
    Array[Array[String]] samples_statuses_and_vcfs
    String output_tsv = "status.tsv"
    RunEnv runenv
  }

  command <<<
    mv write_tsv(~{samples_statuses_and_vcfs}) ~{output_tsv}
  >>>

  output {
    File output_tsv = glob(output_tsv)[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks : runenv.disks
  }
}
