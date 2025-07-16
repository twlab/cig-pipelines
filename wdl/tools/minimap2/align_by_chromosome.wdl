version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/minimap2/align_chromosome.wdl"
import "wdl/tasks/samtools/fai.wdl"
import "wdl/tasks/misc/tar.wdl"

workflow mm2_align_by_chromosome {
  meta {
    author: "Eddie Belter"
    version: "0.1"
    description: "Align Individual Chromosome FASTAs of Query vs Reference"
  }

  input {
    File query_fasta
    File query_fai
    File ref_fasta
    File ref_fai
    String mm2_params
    String mm2_docker
    Int mm2_cpu
    Int mm2_memory
  }

  RunEnv mm2_runenv = {
    "docker": mm2_docker,
    "cpu": mm2_cpu,
    "memory": mm2_memory,
    "disks": 20,
  }

  # Get the Chromosomes from the FAI
  call fai.get_chromosomes { input:
    fai=query_fai,
    runenv=mm2_runenv,
  }

  # Align Each QUERY CHR FASTA to Its REF Mate
  scatter ( chr in get_chromosomes.chromosomes ) {
    call align_chromosome.run_align_chromosome as align { input:
      chr=chr,
      query_fasta=query_fasta,
      query_fai=query_fai,
			ref_fasta=ref_fasta,
      ref_fai=ref_fai,
      params=mm2_params,
      runenv=mm2_runenv
    }
  }

  # TAR the Alignments
  String query_name = basename(query_fasta, ".fasta")
  String ref_name = basename(ref_fasta, ".fasta")
  call tar.run_tar as tarred_alignments { input:
    name="~{ref_name}.~{query_name}.chr-mm2-aln",
    files=align.paf,
    runenv=mm2_runenv,
  }

  output {
    File alignments = tarred_alignments.tar_file
  }
}
