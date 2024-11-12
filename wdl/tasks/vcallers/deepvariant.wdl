version development

import "../../structs/runenv.wdl"

task run_deepvariant {
  input {
    String sample
    File bam
    File bai
    Directory reference_path
    RunEnv runenv
  }

  String output_vcf = "${sample}.vcf.gz"
  String output_gvcf = "${sample}.g.vcf.gz"
  Int dv_cpu = runenv.cpu - 1
  command <<<
    set -ex
    ln ~{bam} ~{basename(bam)}
    ln ~{bai} ~{basename(bai)}
    reference_fasta=$(find ~{reference_path} -name \*.fasta)
    customized_model_param=""
    model_name=$(find ~{reference_path} -name model\* | xargs -I% basename % | awk -F. '{print $1}' | sort -u | head -1)
    if [[ ! -z "${model_name}" ]]; then
      customized_model_param="--customized_model=~{reference_path}/${model_name}"
    fi
    /opt/deepvariant/bin/run_deepvariant \
      --model_type=WGS \
      --ref=${reference_fasta} \
      ${customized_model_param} \
      --reads=~{basename(bam)} \
      --output_vcf=~{output_vcf} \
      --output_gvcf=~{output_gvcf} \
      --num_shards=~{dv_cpu}
    set +e
    printf "Validating VCF: %s\n" "~{output_vcf}" 1>&2
    bcftools view "~{output_vcf}" > /dev/null
    rv=$?
    test "${rv}" != "0" && ( printf "VCF is corrupted, exiting.\n" 1>&2; exit "${rv}" )
    printf "VCF PASS\n" 1>&2
    printf "Validating GVCF: %s\n" "~{output_gvcf}" 1>&2
    bcftools view "~{output_gvcf}" > /dev/null
    rv=$?
    test "${rv}" != "0" && ( printf "GVCF is corrupted, exiting.\n" 1>&2; exit "${rv}" )
    printf "GVCF PASS\n" 1>&2
  >>>

  output {
    String vcf = glob("~{output_vcf}")[0]
    String vcf_tbi = glob("~{output_vcf}.tbi")[0]
    String gvcf = glob("~{output_gvcf}")[0]
    String gvcf_tbi = glob("~{output_gvcf}.tbi")[0]
    String report = glob("*.html")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks : runenv.disks
  }
}

task run_dv {
  input {
    String sample
    File bam
    File bai
    File reference_fasta
    File reference_fai
    String model_type
    RunEnv runenv
  }

  String output_vcf = "${sample}.vcf.gz"
  Int dv_cpu = runenv.cpu - 1
  # MALE
  # --haploid_contigs=chrX,chrY,chrY_KI270740v1_random
  command <<<
    ln ~{bam} ~{basename(bam)}
    ln ~{bai} ~{basename(bai)}
    /opt/deepvariant/bin/run_deepvariant \
      --model_type=~{model_type} \
      --ref=~{reference_fasta} \
      --reads=~{basename(bam)} \
      --output_vcf=~{sample}.vcf.gz \
      --output_gvcf=~{sample}.g.vcf.gz \
      --num_shards=~{dv_cpu}
  >>>

  output {
    File vcf = glob("~{sample}.vcf.gz")[0]
    File vcf_tbi = glob("~{sample}.vcf.gz.tbi")[0]
    File gvcf = glob("~{sample}.g.vcf.gz")[0]
    File gvcf_tbi = glob("~{sample}.g.vcf.gz.tbi")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks : runenv.disks
  }
}
