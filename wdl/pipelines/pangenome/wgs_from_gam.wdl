version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/gatk/realigner_target_creator.wdl"
import "wdl/tasks/pangenome/extract_ref.wdl"
import "wdl/tasks/realign/abra2.wdl"
import "wdl/tasks/realign/freebayes.wdl"
import "wdl/tasks/samtools.wdl"
import "wdl/tasks/vcallers/deepvariant.wdl"
import "wdl/tasks/vg/paths.wdl"
import "wdl/tasks/vg/stats.wdl"
import "wdl/tasks/vg/surject.wdl"

workflow pangenome_wgs_from_gam {
  meta {
    author: "Eddie Belter"
    version: "1.2"
    description: "Run the WGS pipeline astarting with a gam from vg giraffe to surject alignments bam, left/realign indels, and call variants with deepvariant."
  }

  input {
    String sample
    String reference_name
    File gam
    File gbz
    Int targets_expansion_bases = 160
    # dockers
    String abra2_docker = "mgibio/abra2:v2.24-focal"
    String deepvariant_docker = "google/deepvariant:1.5.0"
    String freebayes_docker = "mgibio/freebayes:1.3.6-focal"
    String gatk_docker = "broadinstitute/gatk3@sha256:5ecb139965b86daa9aa85bc531937415d9e98fa8a6b331cb2b05168ac29bc76b" #"broadinstitute/gatk:4.3.0.0"
    String samtools_docker = "mgibio/samtools:1.15.1-buster"
    String vg_docker = "quay.io/vgteam/vg:v1.48.0" #"quay.io/vgteam/vg@sha256:62a1177ab6feb76de6a19af7ad34352bea02cab8aa2996470d9d2b40b3190fe8"
  }

  # RunEnvs in order of usage
  RunEnv runenv_vg = {
    "docker": vg_docker,
    "cpu": 8,
    "memory": 64,
    "disks": 20,
  }

  RunEnv samtools_runenv = {
    "docker": samtools_docker,
    "cpu": 1,
    "memory": 4,
    "disks": 20,
  }

  RunEnv freebayes_renenv = {
    "docker": freebayes_docker,
    "cpu": 1,
    "memory": 4,
    "disks": 20,
  }

  RunEnv gatk_renenv = {
    "docker": gatk_docker,
    "cpu": 4,
    "memory": 24,
    "disks": 20,
  }

  RunEnv abra2_renenv = {
    "docker": abra2_docker,
    "cpu": 2,
    "memory": 20,
    "disks": 20,
  }

  RunEnv dv_runenv = {
    "docker": deepvariant_docker,
    "cpu": 20,
    "memory": 96,
    "disks": 20,
  }

  call stats.run_stats as vg_stats { input:
    gam=gam,
    runenv=runenv_vg,
  }

  call extract_ref.run_extract_ref as reference { input:
    name=reference_name,
    gbz=gbz,
    runenv=runenv_vg,
  }

  call surject.run_surject { input:
    gam=gam,
    sample=sample,
    library=sample+"-lib1",
    gbz=gbz,
    paths_list=reference.paths_list,
    runenv=runenv_vg,
  }

  call samtools.sort as samtools_sort { input:
    bam=run_surject.bam,
    runenv=samtools_runenv,
  } 

  call freebayes.run_left_align_bam as left_align { input:
    in_bam_file=samtools_sort.sorted_bam,
    in_reference_file=reference.fasta,
    in_reference_index_file=reference.fai,
    runenv=freebayes_renenv,
  } 

  call samtools.index as left_align_index { input:
    bam=left_align.output_bam_file,
    runenv=samtools_runenv,
  } 

  call realigner_target_creator.run_realigner_target_creator as target_creator { input:
    bam=left_align.output_bam_file,
    bai=left_align_index.bai,
    expand_bases=targets_expansion_bases,
    reference_fasta=reference.fasta,
    reference_fai=reference.fai,
    reference_dict=reference.dict,
    runenv=gatk_renenv,
  } 

  call abra2.run_realigner as realign { input:
    in_bam_file=left_align.output_bam_file,
    in_bam_index_file=left_align_index.bai,
    in_target_bed_file=target_creator.expanded_targets,
    in_reference_file=reference.fasta,
    in_reference_index_file=reference.fai,
    runenv=abra2_renenv,
  } 

  call deepvariant.run_deepvariant as dv { input:
    sample=sample,
    bam=realign.indel_realigned_bam,
    bai=realign.indel_realigned_bam_index,
    ref_fasta=reference.fasta,
    ref_fai=reference.fai,
    ref_dict=reference.dict,
    runenv=dv_runenv,
  }

  call samtools.stat as samtools_stat { input:
    bam=realign.indel_realigned_bam,
    runenv=samtools_runenv,
  } 

  output {
    File gam_stats = vg_stats.stats
    File bam = realign.indel_realigned_bam
    File bai = realign.indel_realigned_bam_index
    File bam_stats = samtools_stat.stats
    File dv_vcf = dv.vcf
  }
}
