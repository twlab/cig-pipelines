# Giraffe Pipeline

Map WGS data to the pangenome graph, and call variants with deep variant. Produces a GAM, sorted BAM & BAI, plus VCF & stats.

## Pipeline Files
* giraffe.wdl          - WDL pipeline
* giraffe.inputs.json  - pipeline inputs with place holders
* giraffe.outputs.yaml - steps and outputs to be copied after pipeline run
* giraffe.imports.zip  - imports used in the WDL
* giraffe.doc.md       - this file, documenting the pipeline

## Inputs
* sample [String] - sample name for outputs
* fastqs [File] - an array of read1 and read2 fastqs
* gbz [File] - giraffe pangenome GBZ

## Steps
### Untar the BWA Reference
#### input
* idx [inputs.idx]
#### output
* reference - untarred BWA idx

### Align with BWA MEM
Align sets of read 1 & 2 fastqs
#### input
* name [workflow inputs]
* fastqs [workflow inputs]
* reference [output from untar idx]
####output:
* bam

### Merge Bams
#### input
* name [workflow inputs]
* bams [output from align]
####output:
* merged_bam

### Collect Samtools Stat
#### input
* bam [output from align]
#### output
* stats [samtools stats file]

### Sort the BAM by Coordinates
The bam needs to be sorted by coordinate to call variants
#### input
* bam [output from merge]
#### output
* sorted_bam

### Picard Mark Duplicates
#### input
* name [workflow inputs]
* bam [output bam from samtools sort]
#### output
* dedup_bam
* metrics

### BQSR
[GATK BQSR doc](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-)
#### input
* bam [output bam from picard markdup]
* reference [output from untar idx]
* known_sites [workflow input]
#### output
* recalibrated bam
 
### Haplotype Caller for SNPs and INDELs
[GATK Haplotype Caller doc](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller)
#### input
* bam [otuput bam from gatk bqsr]
* reference [output from untar idx]
#### output
* vcf

## Outputs
* bam [from gatk_bqsr.recal_bam]
* vcf [from haplotype caller]
* dedup_metrics [from markdup]
* stats [from stats]
version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/samtools.wdl"
import "wdl/tasks/vg/giraffe.wdl"
import "wdl/tasks/vg/stats.wdl"
import "wdl/tasks/vg/surject.wdl"

workflow giraffe_pipeline {

  input {
    String name
    Array[File] fastqs 
    File min
    File dist
    File gbz
    String docker = "quay.io/vgteam/vg:v1.48.0" #"quay.io/vgteam/vg@sha256:62a1177ab6feb76de6a19af7ad34352bea02cab8aa2996470d9d2b40b3190fe8"
    Int cpu = 32
    Int memory = 500
  }

  RunEnv giraffe_runenv = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
    "disks": 20,
  }

  call giraffe.run_giraffe { input:
    name=name,
    fastqs=fastqs,
    min=min,
    dist=dist,
    gbz=gbz,
    runenv=giraffe_runenv,
  }

  RunEnv runenv_vg = {
    "docker": docker,
    "cpu": 4,
    "memory": 20,
    "disks": 20,
  }

  call stats.run_stats { input:
    gam=run_giraffe.gam,
    runenv=runenv_vg,
  }

  call surject.run_surject { input:
    gam=run_giraffe.gam,
    gbz=gbz,
    runenv=runenv_vg,
  }

  call samtools.stat as samtools_stat { input:
    bam=run_surject.bam,
    runenv=runenv_vg, # vg 1.48.0 docker has samtools 1.10
  }

  output {
    File gam = run_giraffe.gam
    File gam_stats = run_stats.stats
    File bam = run_surject.bam
    File bam_stats = samtools_stat.stats
  }
}
