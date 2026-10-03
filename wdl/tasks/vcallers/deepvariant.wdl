version development

import "../../structs/runenv.wdl"

task run_deepvariant {
  input {
    String sample
    File bam
    File bai
    File ref_fasta
    File ref_fai
    File ref_dict
    Boolean generate_gvcf = true
    String model_type = "WGS"
    RunEnv runenv
  }

  String output_vcf = "~{sample}.vcf.gz"
  String output_gvcf = "~{sample}.g.vcf.gz"
  command <<<
    set -ex
    ln ~{bam} ~{basename(bam)}
    ln ~{bai} ~{basename(bai)}
    ln ~{ref_fasta} ~{basename(ref_fasta)}
    ln ~{ref_fai} ~{basename(ref_fai)}
    ln ~{ref_dict} ~{basename(ref_dict)}

    /opt/deepvariant/bin/run_deepvariant \
      --sample_name=~{sample} \
      --model_type=~{model_type} \
      --ref=~{basename(ref_fasta)} \
      --reads=~{basename(bam)} \
      --output_vcf=~{output_vcf} \
      ~{if (generate_gvcf) then "--output_gvcf=~{output_gvcf}}" else ""} \
      --num_shards=~{runenv.cpu}
    set +e
    printf "Validating VCF: %s\n" "~{output_vcf}" 1>&2
    bcftools view "~{output_vcf}" > /dev/null
    rv=$?
    test "${rv}" != "0" && ( printf "VCF is corrupted, exiting.\n" 1>&2; exit "${rv}" )
    printf "VCF PASS\n" 1>&2
    if test -e "~{output_gvcf}"; then
      printf "Validating GVCF: %s\n" "~{output_gvcf}" 1>&2
      bcftools view "~{output_gvcf}" > /dev/null
      rv=$?
      test "${rv}" != "0" && ( printf "GVCF is corrupted, exiting.\n" 1>&2; exit "${rv}" )
      printf "GVCF PASS\n" 1>&2
    fi
  >>>

  output {
    File vcf = output_vcf
    File vcf_tbi = "~{output_vcf}.tbi"
    File? gvcf = output_gvcf
    File? gvcf_tbi = "~{output_gvcf}.tbi"
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks : runenv.disks
  }
}
