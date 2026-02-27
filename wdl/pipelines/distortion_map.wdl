version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/bwa/align.wdl"
import "wdl/tasks/bwa/idx.wdl"
import "wdl/tasks/bed/bedtools.wdl"
import "wdl/tasks/distortion_map/db.wdl"
import "wdl/tasks/distortion_map/count_matrices.wdl"
import "wdl/tasks/distortion_map/coverage.wdl"
import "wdl/tasks/distortion_map/metrics.wdl"
import "wdl/tasks/distortion_map/intervals.wdl"
import "wdl/tasks/distortion_map/wgsim.wdl"
import "wdl/tasks/misc/cat.wdl"
import "wdl/tasks/minimap2/align.wdl" as minimap2_align
import "wdl/tasks/minimap2/liftover.wdl"
import "wdl/tasks/samtools/faidx.wdl" as samtools_faidx

workflow distortion_map {
  input {
    File target_fasta           # target (ne: reference) full FASTA
    File target_fasta_gzi       # target GZI
    File query_fasta            # query full FASTA
    File query_fasta_gzi        # query GZI
    Array[Int] wgsim_coverages  # List coverages to run, must be different to avoid cromwell cache hit
    Float wgsim_base_error
    Int wgsim_out_distance
    Int wgsim_stdev
    Int wgsim_read1_length
    Int wgsim_read2_length
    Int interval_window_length
    Int interval_window_stride
    Float wgsim_mutation_rate
    Float wgsim_fraction_indels
    Float wgsim_prob_indel_extentsion
    Int wgsim_seed
    # Resources
    String bwa_docker
    Int bwa_cpu
    Int bwa_memory
    String dm_docker
    Int dm_generate_count_mtx_cpu
    Int dm_generate_count_mtx_memory
    String minimap2_docker
    String wgsim_docker
    Int wgsim_cpu
    Int wgsim_memory
  }

  # RunEnvs
  # DM Run Envs
  # *can thread
  # 
  # Prepare Target & Query
  #  extract chr    dm 1 cpu / 4 G
  #  bwa idx        bwa ?
  #  minimap2*      dm 8 cpu / 64 G
  # Gen Count Mtx
  #  read cov    
  #  simulate
  #  read pos
  #  bedtools       dm 1 cpu / 4 G
  #  liftover       dm 1 cpu / 48 G
  #  load db        dm 1 cpu / 48 G
  #  intervals      em 1 cpu / 24 G
  #  generate mtx*  16 cpu / 64 G
  # Merge, Normalize, & Gen Metrics
  #  merge mtx      dm 1 cpu / 24 G
  #  normalize mtx  dm 1 cpu / 24 G
  #  calculate mtx  dm 1 cpu / 24 G

  # Reviewed Runenvs
  RunEnv bwa_runenv = {
    "docker": bwa_docker,
    "cpu": bwa_cpu,
    "memory": bwa_memory,
    "disks": 20,
  }

  RunEnv dm_runenv_1cpu_4G = {
    "docker": dm_docker,
    "cpu": 1,
    "memory": 4,
    "disks": 20,
  }

  RunEnv dm_runenv_1cpu_24G = {
    "docker": dm_docker,
    "cpu": 1,
    "memory": 24,
    "disks": 20,
  }

  RunEnv dm_runenv_1cpu_48G = {
    "docker": dm_docker,
    "cpu": 1,
    "memory": 48,
    "disks": 20,
  }

  RunEnv dm_runenv_1cpu_72G = {
    "docker": dm_docker,
    "cpu": 1,
    "memory": 72,
    "disks": 20,
  }

  RunEnv dm_runenv_minimap2 = {
    "docker": dm_docker,
    "cpu": 6,
    "memory": 48,
    "disks": 20,
  }

  RunEnv dm_runenv_generate_count_mtx = {
    "docker": dm_docker,
    "cpu": dm_generate_count_mtx_cpu,
    "memory": dm_generate_count_mtx_memory,
    "disks": 20,
  }

  RunEnv minimap2_runenv_1cpu_48G = {
    "docker": minimap2_docker,
    "cpu": 1,
    "memory": 48,
    "disks": 20,
  }

  RunEnv wgsim_runenv = {
    "docker": wgsim_docker,
    "cpu": wgsim_cpu,
    "memory": wgsim_memory,
    "disks": 20,
  }

  # Get Chromosomes from Query FASTA
  call samtools_faidx.extract_chromosome_names as chromosomes { input:
    fasta=query_fasta,
    runenv=dm_runenv_1cpu_4G,
  }

  # Prepare Query and Target Resources
  #  * BWA Indexes for Query Chromosomes
  #  * BWA Indexes for Target Chromosomes
  #  * Minimap2 PAF for Query vs. Target
  String query_name = basename(query_fasta, ".fasta.gz")
  String target_name = basename(target_fasta, ".fasta.gz")
  scatter (chr in chromosomes.names) {
    # Query - extract chr fasta & build index
    call samtools_faidx.extract_chromosome as query_chr_fasta { input:
      fasta=query_fasta,
      fai=query_fasta_gzi,
      chr=chr,
      runenv=dm_runenv_1cpu_4G,
    }
    call idx.run_build_idx as query_chr { input:
      name="~{query_name}_~{chr}",
      fasta=query_chr_fasta.chr_fasta,
      runenv=bwa_runenv,
    }

    # Target - extract chr fasta & build index
    call samtools_faidx.extract_chromosome as target_chr_fasta { input:
      fasta=target_fasta,
      fai=target_fasta_gzi,
      chr=chr,
      runenv=dm_runenv_1cpu_4G,
    }
    call idx.run_build_idx as target_chr { input:
      name="~{target_name}_~{chr}",
      fasta=target_chr_fasta.chr_fasta,
      runenv=bwa_runenv,
    }

    # MINIMAP2 - Query to Target PAF
    call minimap2_align.run_align as query_to_ref_paf { input:
      query_fasta=query_chr_fasta.chr_fasta,
      target_fasta=target_chr_fasta.chr_fasta,
      output_fn="~{target_name}.~{query_name}.paf",
      params="-x asm5 -L -c --cs=long",
      runenv=dm_runenv_minimap2,
    }
  }

  # Combine the Target Sizes and Create Intervals
  call cat.cat as target_sizes { input:
    files=target_chr.sizes,
    out="~{target_name}.sizes",
    runenv=dm_runenv_1cpu_4G,
  }

  # Run Each Chromosome Through Each Coverage
  Array[Int] chromosomes_i = range(length(chromosomes.names))
  scatter (i in chromosomes_i) {
    String chr = chromosomes.names[i]
    scatter (wgsim_coverage in wgsim_coverages) {
      # Simulate reads from query chromosome
      call wgsim.calc_read_pairs_needed { input:
        fasta=query_chr_fasta.chr_fasta[i],
        coverage=wgsim_coverage,
        read_length=wgsim_read1_length,
        runenv=dm_runenv_1cpu_4G,
      }
      call wgsim.run_wgsim { input:
        fasta=query_chr_fasta.chr_fasta[i],
        base_error=wgsim_base_error,
        out_distance=wgsim_out_distance,
        stdev=wgsim_stdev,
        number_pairs=calc_read_pairs_needed.read_pairs_needed,
        read1_length=wgsim_read1_length,
        read2_length=wgsim_read2_length,
        mutation_rate=wgsim_mutation_rate,
        fraction_indels=wgsim_fraction_indels,
        prob_indel_extentsion=wgsim_prob_indel_extentsion,
        seed=wgsim_seed+wgsim_coverage,
        runenv=wgsim_runenv,
      }
  
      # Extract Source Locations from FASTQs
      call wgsim.extract_source_positions { input:
        fastqs=[run_wgsim.simulated_r1_fastq, run_wgsim.simulated_r2_fastq],
        runenv=wgsim_runenv,
      }
      call liftover.run_liftover as liftover_source_positions_to_ref { input:
        paf=query_to_ref_paf.paf[i],
        bed=extract_source_positions.source_positions,
        mapping_qual=5,
        alignment_length=50000,
        max_seq_divergence=1,
        runenv=minimap2_runenv_1cpu_48G,
      }
  
      # Align Sim Reads to REF then Convert Alignments to BED
      String sample = "~{query_name}_~{chr}"
      call align.run_bwa_mem as align_to_ref { input:
        sample=sample,
        library=sample+"-lib1",
        rg_id=sample+"-lib1",
        platform_unit="ILLUMINA",
        fastqs=[run_wgsim.simulated_r1_fastq, run_wgsim.simulated_r2_fastq],
        idx_files=[target_chr.FASTA[i], target_chr.amb[i], target_chr.ann[i], target_chr.bwt[i], target_chr.pac[i], target_chr.sa[i]], 
        runenv=bwa_runenv,
      }
      call bedtools.run_bam_to_bed as bam2bed_ref { input:
        bam=align_to_ref.bam,
        runenv=dm_runenv_1cpu_24G,
      }
  
      # Align Sim Reads to QUERY then Convert Alignments to BED & Lift Over
      call align.run_bwa_mem as align_to_query { input:
        sample=sample,
        library=sample+"-lib1",
        rg_id=sample+"-lib1",
        platform_unit="ILLUMINA",
        fastqs=[run_wgsim.simulated_r1_fastq, run_wgsim.simulated_r2_fastq],
        idx_files=[query_chr.FASTA[i], query_chr.amb[i], query_chr.ann[i], query_chr.bwt[i], query_chr.pac[i], query_chr.sa[i]],
        runenv=bwa_runenv,
      }
      call bedtools.run_bam_to_bed as bam2bed_query { input:
        bam=align_to_query.bam,
        runenv=dm_runenv_1cpu_24G,
      }
      call liftover.run_liftover as liftover_query_alignments_to_ref { input:
        paf=query_to_ref_paf.paf[i],
        bed=bam2bed_query.bedfile,
        mapping_qual=5,
        alignment_length=50000,
        max_seq_divergence=1,
        runenv=minimap2_runenv_1cpu_48G,
      }

      call samtools_faidx.extract_chromosome_size as chr_size { input:
        chr=chr,
        fai=query_chr.fai[i],
        runenv=dm_runenv_1cpu_4G,
      }
      RunEnv dm_load_db_runenv = if ( chr_size.size >= 175000000 ) then dm_runenv_1cpu_72G else dm_runenv_1cpu_48G

      call db.load_db { input:
        source_positions=extract_source_positions.source_positions,
        lifted_source=liftover_source_positions_to_ref.bedfile,
        aligned_ref=bam2bed_ref.bedfile,
        aligned_source=bam2bed_query.bedfile,
        lifted_aligned_source=liftover_query_alignments_to_ref.bedfile,
        runenv=dm_load_db_runenv,
      }
  
      call intervals.create_intervals as intervals { input:
        db=load_db.db,
        reference_sizes=target_sizes.concatenated_file,
        window_length=interval_window_length,
        window_stride=interval_window_stride,
        runenv=dm_runenv_1cpu_24G,
      }
  
      call count_matrices.generate as count_mtx { input:
        db=load_db.db,
        reference_intervals=intervals.reference_intervals,
        runenv=dm_runenv_generate_count_mtx,
      }
    }
  }

  call count_matrices.merge as merge_ref_matrices { input:
    matrices=flatten(count_mtx.aligned_reference_count_matrix),
    output_fn="merged.ref.mtx",
    runenv=dm_runenv_1cpu_24G,
  }

  call count_matrices.merge as merge_src_matrices { input:
    matrices=flatten(count_mtx.lifted_aligned_source_count_matrix),
    output_fn="merged.src.mtx",
    runenv=dm_runenv_1cpu_24G,
  }

  call count_matrices.normalize as normalize_merged_ref_matrix { input:
    matrix=merge_ref_matrices.matrix,
    runenv=dm_runenv_1cpu_24G,
  }

  call count_matrices.normalize as normalize_merged_src_matrix { input:
    matrix=merge_src_matrices.matrix,
    runenv=dm_runenv_1cpu_24G,
  }

  call metrics.calculate as calculate_metrics { input:
    normalized_aligned_reference_matrix=normalize_merged_ref_matrix.normalized_matrix,
    normalized_lifted_aligned_source_matrix=normalize_merged_src_matrix.normalized_matrix,
    interval_mapping=select_first(flatten(count_mtx.interval_mapping)),
    runenv=dm_runenv_1cpu_24G,
  }
}
