version development

import "wdl/structs/runenv.wdl"

workflow minimap2 {
  input {
    File fastq
    Directory reference
		File junctions
    String output_prefix
    String read_type
    String docker
    Int cpu
    Int memory
  }

  RunEnv runenv = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
    "disks": 20,
  }

  call run_align as minimap_align { input:
    fastq=fastq,
    reference=reference,
    junctions=junctions,
    output_prefix=output_prefix,
    read_type=read_type,
    runenv=runenv,
  }
}
 
task run_align {
  input {
    File fastq
    String read_type
    Directory reference 
    File junctions
    String output_prefix 
    RunEnv runenv
  }

  String output_bam = "~{output_prefix}.bam"
  String stat_fn = "~{output_bam}.stat"
  # FIXME ENCODE uses different params for PacBio/Ont
  # pacbio: -ax splice -uf --secondary=no -C5
  # ont:    -ax splice -uf -k14 \
  command <<<
    set -ex
    reference_fasta=$(find ~{reference} -name \*.fa | head -1)
    minimap2 -t ~{runenv.cpu} -ax splice -k14 -uf --junc-bed ~{junctions} $reference_fasta ~{fastq} | samtools view -q 40 -F2304 -bS - | samtools sort - -O bam -@ ~{runenv.cpu} > ~{output_bam}
    samtools stats ~{output_bam} > ~{stat_fn}
  >>>

  output {
    File output_bam = glob("~{output_bam}")[0]
    File stats = glob("~{stat_fn}")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    #disks : select_first([runenv.disks,"local-disk 100 SSD"])
  }
}
