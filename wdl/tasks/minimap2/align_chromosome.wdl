version development

import "../../structs/runenv.wdl"

task run_align_chromosome {
  input {
    String chr
    File query_fasta
    File query_fai
    File ref_fasta
    File ref_fai
    String params
    RunEnv runenv
  }

  String query_name = basename(query_fasta, ".fasta")
  String ref_name = basename(ref_fasta, ".fasta")
  command <<<
    samtools faidx ~{query_fasta} --fai-idx ~{query_fai} ~{chr} -o ~{query_name}.~{chr}.fasta
    samtools faidx ~{ref_fasta} --fai-idx ~{ref_fai} ~{chr} -o ~{ref_name}.~{chr}.fasta
    minimap2 ~{params} -o ~{ref_name}.~{query_name}.paf ~{ref_fasta} ~{query_fasta}
  >>>

  output {
    File paf = glob("*.paf")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    #disks : select_first([runenv.disks,"local-disk 100 SSD"])
  }
}
