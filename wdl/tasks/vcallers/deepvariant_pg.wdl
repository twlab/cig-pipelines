version development

import "../../structs/runenv.wdl"

task run_deepvariant_pg {
  input {
    String sample
    File bam
    File bai
    File ref_fasta
    File ref_fai
    File ref_dict
    File gbz
    String model_type = "WGS"
    RunEnv runenv
  }

  String output_vcf = "~{sample}.vcf.gz"
  String output_gvcf = "~{sample}.g.vcf.gz"
  # --gbz_shared_memory_size_gb: Optional. Size of the shared memory region to create. (default: '12') (an integer)
  command <<<
    set -ex
    ln ~{bam} ~{basename(bam)}
    ln ~{bai} ~{basename(bai)}
    ln ~{ref_fasta} ~{basename(ref_fasta)}
    ln ~{ref_fai} ~{basename(ref_fai)}
    ln ~{ref_dict} ~{basename(ref_dict)}

    /opt/deepvariant/bin/run_pangenome_aware_deepvariant \
      --model_type=~{model_type} \
      --ref=~{basename(ref_fasta)} \
      --pangenome ~{gbz} \
      --reads=~{basename(bam)} \
      --sample_name_reads=~{sample} \
      --output_vcf=~{output_vcf} \
      --output_gvcf=~{output_gvcf} \
      --num_shards=~{runenv.cpu}

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
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks: runenv.disks
    shm_size: "12g" # docker shared memory size
  }
}
