version development

import "../../structs/runenv.wdl"

task run_verifybamid {
  input {
    File input_bam
    RunEnv runenv
  }

  String output_file = sub(basename(input_vcf), ".vcf.gz$", ".checksex.tsv")
  String input_vcf_bn = basename(input_vcf)
  String input_vcf_tbi_bn = basename(input_vcf_tbi)
  command <<<
    set -x

verifybamid2 --SVDPrefix /opt/conda/share/verifybamid2-2.0.1-12/resource/1000g.phase3.100k.b38.vcf.gz.dat --Reference /storage2/fs1/epigenome/Active/ctomlins/VerifyBAMID2_Reference_Files/GRCh38_full_analysis_set_plus_decoy_hla.fasta --BamFile /storage2/fs1/epigenome/Active/shared_smaht/Human_Contamination_Evaluation_SR004399/Alignments_GRCh38/LIB035789-DIL01_227KVNLT4_S3_L002.GRCh38.sorted.bam

    ln -s ~{input_vcf} .
    ln -s ~{input_vcf_tbi} .
    smaht tools checksex ~{input_vcf_bn} ~{output_file} -p ~{seq_platform}
  >>>

  output {
    File output_file = glob("~{output_file}")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "${runenv.memory} GB"
    #disks: runenv.disks
  }
}
