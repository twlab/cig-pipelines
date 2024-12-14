version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/bwa/align.wdl"
import "wdl/tasks/bwa/idx.wdl"
import "wdl/tasks/bed/bedtools.wdl"
import "wdl/tasks/distortion_map/db.wdl"
import "wdl/tasks/distortion_map/count_matrices.wdl"
import "wdl/tasks/distortion_map/coverage.wdl"
import "wdl/tasks/distortion_map/intervals.wdl"
import "wdl/tasks/distortion_map/wgsim.wdl"
import "wdl/tasks/minimap2/liftover.wdl"
import "wdl/tasks/samtools/split.wdl"

workflow distortion_map {
    input {
      String sample
      File query_idx          # tar file with ref fasta, fai, and aligner index
      File reference_idx      # tar file with ref fasta, fai, and aligner index
      File query_to_ref_paf
      Array[String] chrs
      Float wgsim_base_error
      Int wgsim_out_distance
      Int wgsim_stdev
      Int wgsim_number_pairs
      Int wgsim_read1_length
      Int wgsim_read2_length
      Float wgsim_mutation_rate
      Float wgsim_fraction_indels
      Float wgsim_prob_indel_extentsion
      Int wgsim_seed
      # Resources
      String bwa_docker
      Int bwa_cpu
      Int bwa_memory
      String bedtools_docker
      Int bedtools_cpu
      Int bedtools_memory
      String distortion_map_docker
      Int distortion_map_cpu
      Int distortion_map_memory
      String minimap2_docker
      Int minimap2_cpu
      Int minimap2_memory
      String samtools_docker
      Int samtools_cpu
      Int samtools_memory
      String wgsim_docker
      Int wgsim_cpu
      Int wgsim_memory
      String utils_docker
      Int utils_cpu
      Int utils_memory
    }

  # RunEnvs in order of usage
  RunEnv utils_runenv = {
    "docker": utils_docker,
    "cpu": utils_cpu,
    "memory": utils_memory,
    "disks": 20,
  }

  RunEnv samtools_runenv = {
    "docker": samtools_docker,
    "cpu": samtools_cpu,
    "memory": samtools_memory,
    "disks": 20,
  }

  RunEnv wgsim_runenv = {
    "docker": wgsim_docker,
    "cpu": wgsim_cpu,
    "memory": wgsim_memory,
    "disks": 20,
  }

  RunEnv bwa_runenv = {
    "docker": bwa_docker,
    "cpu": bwa_cpu,
    "memory": bwa_memory,
    "disks": 20,
  }

  RunEnv bedtools_runenv = {
    "docker": bedtools_docker,
    "cpu": bedtools_cpu,
    "memory": bedtools_memory,
    "disks": 20,
  }

  RunEnv minimap2_runenv = {
    "docker": minimap2_docker,
    "cpu": minimap2_cpu,
    "memory": minimap2_memory,
    "disks": 20,
  }

  RunEnv distortion_map_runenv = {
    "docker": distortion_map_docker,
    "cpu": distortion_map_cpu,
    "memory": distortion_map_memory,
    "disks": 20,
  }

  RunEnv distortion_map_threaded_runenv = {
    "docker": distortion_map_docker,
    "cpu": distortion_map_cpu * 4,
    "memory": distortion_map_memory * 2,
    "disks": 20,
  }

  # Untar the Indexes
  call idx.run_untar_idx as reference { input:
    idx=reference_idx,
    runenv=utils_runenv,
  }

  call idx.run_untar_idx as query { input:
    idx=query_idx,
    runenv=utils_runenv,
  }

  # Split the QUERY FASTA by Chromosome
  call split.run_split_by_chromosome as splitter { input:
    fasta=query.fasta,
    fai=query.fai,
    chrs=chrs,
    runenv=samtools_runenv,
  }

  # Process Each QUERY Chromosome
  scatter (chromosome_fasta in splitter.chromosome_fastas) {
    # Simulate reads from query chromosome
    call wgsim.run_wgsim { input:
      fasta=chromosome_fasta,
      base_error=wgsim_base_error,
      out_distance=wgsim_out_distance,
      stdev=wgsim_stdev,
      number_pairs=wgsim_number_pairs,
      read1_length=wgsim_read1_length,
      read2_length=wgsim_read2_length,
      mutation_rate=wgsim_mutation_rate,
      fraction_indels=wgsim_fraction_indels,
      prob_indel_extentsion=wgsim_prob_indel_extentsion,
      seed=wgsim_seed,
      runenv=wgsim_runenv,
    }

    # Extract Source Locations from FASTQs
    call wgsim.extract_source_positions { input:
      fastqs=[run_wgsim.simulated_r1_fastq, run_wgsim.simulated_r2_fastq],
      runenv=wgsim_runenv,
    }
    call liftover.run_liftover as liftover_source_positions_to_ref { input:
      paf=query_to_ref_paf,
      bed=extract_source_positions.source_positions,
      mapping_qual=5,
      alignment_length=50000,
      max_seq_divergence=1,
      runenv=minimap2_runenv,
    }

    # Align Sim Reads to REF
    call align.run_bwa_mem as align_to_ref { input:
      sample=sample,
      library=sample+"-lib1",
      fastqs=[run_wgsim.simulated_r1_fastq, run_wgsim.simulated_r2_fastq],
      reference=reference.path,
      runenv=bwa_runenv,
    }
    call bedtools.run_bam_to_bed as bam2bed_ref { input:
      bam=align_to_ref.bam,
      runenv=bedtools_runenv,
    }

    # Align Sim Reads to QUERY
    call align.run_bwa_mem as align_to_query { input:
      sample=sample,
      library=sample+"-lib1",
      fastqs=[run_wgsim.simulated_r1_fastq, run_wgsim.simulated_r2_fastq],
      reference=query.path,
      runenv=bwa_runenv,
    }
    call bedtools.run_bam_to_bed as bam2bed_query { input:
      bam=align_to_query.bam,
      runenv=bedtools_runenv,
    }
    call liftover.run_liftover as liftover_query_alignments_to_ref { input:
      paf=query_to_ref_paf,
      bed=bam2bed_query.bedfile,
      mapping_qual=5,
      alignment_length=50000,
      max_seq_divergence=1,
      runenv=minimap2_runenv,
    }

    # Load the Database
    call db.load_db { input:
      source_positions=extract_source_positions.source_positions,
      lifted_source=liftover_source_positions_to_ref.bedfile,
      aligned_ref=bam2bed_ref.bedfile,
      aligned_source=bam2bed_query.bedfile,
      lifted_aligned_source=liftover_query_alignments_to_ref.bedfile,
      runenv=distortion_map_runenv,
    }

    # Create the Intervals
    call intervals.create_intervals { input:
      db=load_db.db,
      reference_sizes=reference.sizes,
      window_length=100000,
      window_stride=100000,
      runenv=distortion_map_runenv,
    }
    #File reference_intervals = glob("reference_intervals.tsv")[0]
    #File simulated_intervals = glob("simulated_intervals.tsv")[0]
    #File source_sizes = glob("source_sizes.tsv")[0]

    # Generate Coverages
    call coverage.generate_simulated_coverage { input:
      db=load_db.db,
      simulated_intervals=create_intervals.simulated_intervals,
      batch_size=100000,
      runenv=distortion_map_runenv,
    }

    call coverage.generate_simulated_no_lift_over_coverage { input:
      db=load_db.db,
      simulated_intervals=create_intervals.simulated_intervals,
      runenv=distortion_map_runenv,
    }

    # Count Matrices
    call count_matrices.generate_count_matrices { input:
      db=load_db.db,
      reference_intervals=create_intervals.reference_intervals,
      runenv=distortion_map_threaded_runenv,
    }
  }
}
