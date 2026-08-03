version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/distortion_map/wgsim.wdl"
import "wdl/tasks/minibwa/align.wdl" as minibwa_align
import "wdl/tasks/minibwa/idx.wdl" as minibwa_idx
import "wdl/tasks/minimap2/align.wdl" as minimap2_align
import "wdl/tasks/misc/tar.wdl"
import "wdl/tasks/refract.wdl"
import "wdl/tasks/samtools/faidx.wdl" as samtools_faidx

workflow morgana_simulate_reads {
  input {
    File ref_idx                # minibwa index tar with FASTA, fai, minibwa index files, etc
    File query_fasta            # query FASTA
    Int coverage
    Float wgsim_base_error
    Int wgsim_out_distance
    Int wgsim_stdev
    Int wgsim_read1_length
    Int wgsim_read2_length
    Float wgsim_mutation_rate
    Float wgsim_fraction_indels
    Float wgsim_prob_indel_extentsion
    Int wgsim_seed
    String morgana_docker
  }

  RunEnv minibwa_idx_runenv = {
    "docker": morgana_docker,
    "cpu": 4,
    "memory": 72,
    "disks": 20,
  }

  RunEnv minibwa_align_runenv = {
    "docker": morgana_docker,
    "cpu": 12,
    "memory": 96,
    "disks": 20,
  }

  RunEnv morgana_runenv_1cpu_4G = {
    "docker": morgana_docker,
    "cpu": 1,
    "memory": 4,
    "disks": 20,
  }

  RunEnv minimap2_runenv = {
    "docker": morgana_docker,
    "cpu": 6,
    "memory": 48,
    "disks": 20,
  }

  # Untar REF minibwa index
  call minibwa_idx.run_untar_idx as ref { input:
    idx=ref_idx,
    runenv=morgana_runenv_1cpu_4G,
  }

  # Create QUERY minibwa index
  String query_name = sub(basename(query_fasta), "\\..*", "")
  call minibwa_idx.run_build_idx as query { input:
    name=query_name,
    fasta=query_fasta,
    runenv=minibwa_idx_runenv,
  }

  # QUERY to REF alignment - minimap2 PAF
  String ref_name = sub(basename(ref.fasta), "\\..*", "")
  call minimap2_align.run_align as query_to_ref_paf { input:
    query=query.FASTA,
    target=ref.fasta,
    output_fn="~{ref_name}.~{query_name}.paf",
    params="-x asm5 -L -c --cs=long",
    runenv=minimap2_runenv,
  }

  # QUERY to REF PAF refract index
  call refract.build_idx as refract_rfx { input:
    paf=query_to_ref_paf.alignments,
    runenv=morgana_runenv_1cpu_4G,
  }

  # Simulate reads from QUERY
  # Calculate number of reads to simulate
  call wgsim.calculate_read_pairs_needed { input:
    fai=query.fai,
    coverage=coverage,
    read_length=wgsim_read1_length,
    runenv=morgana_runenv_1cpu_4G,
  }

  # Simulate reads from QUERY
  call wgsim.run_wgsim { input:
    fasta=query.FASTA,
    base_error=wgsim_base_error,
    out_distance=wgsim_out_distance,
    stdev=wgsim_stdev,
    number_pairs=calculate_read_pairs_needed.read_pairs_needed,
    read1_length=wgsim_read1_length,
    read2_length=wgsim_read2_length,
    mutation_rate=wgsim_mutation_rate,
    fraction_indels=wgsim_fraction_indels,
    prob_indel_extentsion=wgsim_prob_indel_extentsion,
    seed=wgsim_seed,
    runenv=morgana_runenv_1cpu_4G,
  }

  # Align simulated reads to full ref and query
  # Align simulated reads to full REF
  String sample = "~{query_name}"
  call minibwa_align.run_minibwa as align_to_ref { input:
    sample=sample,
    library=sample+"-lib1",
    rg_id=sample+"-lib1",
    platform_unit="ILLUMINA",
    fastqs=[run_wgsim.simulated_r1_fastq, run_wgsim.simulated_r2_fastq],
    idx_files=[ref.fasta, ref.l2b, ref.mbw],
    output_fn="~{ref_name}.sorted.bam",
    runenv=minibwa_align_runenv,
  }

  # Align simulated reads to full QUERY
  call minibwa_align.run_minibwa as align_to_query { input:
    sample=sample,
    library=sample+"-lib1",
    rg_id=sample+"-lib1",
    platform_unit="ILLUMINA",
    fastqs=[run_wgsim.simulated_r1_fastq, run_wgsim.simulated_r2_fastq],
    idx_files=[query.FASTA, query.l2b, query.mbw],
    output_fn="~{query_name}.sorted.bam",
    runenv=minibwa_align_runenv,
  }

  # TAR alignments
  call tar.run_tar as tar_simulated_reads { input:
    files=[align_to_ref.output_fn, align_to_query.output_fn, query_to_ref_paf.alignments, refract_rfx.rfx],
    output_file="~{query_name}.~{ref_name}.simreads.tar",
    runenv=morgana_runenv_1cpu_4G,
  }

  output {
    File simulated_reads = tar_simulated_reads.tar_file
  }
}
