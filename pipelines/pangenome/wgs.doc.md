# Pangenome WGS Pipeline

Map WGS data to the pangenome graph using giraffe, and call variants with haplotype caller and deep variant. Produces a GAM, sorted BAM & BAI, plus VCF & stats.

## Pipeline Files
* wgs.wdl          - WDL pipeline
* wgs.inputs.json  - pipeline inputs with place holders
* wgs.outputs.yaml - steps and outputs to be copied after pipeline run
* wgs.imports.zip  - imports used in the WDL
* wgs.doc.md       - this file, documenting the pipeline

## Inputs
* sample [String] - sample name for outputs
* fastqs [File] - an array of read1 and read2 fastqs
* gbz [File] - giraffe pangenome GBZ
* dist - pangenome dist
* min - pangenome mxoin
* reference - linear reference path with FAI and DICT

## Steps
### Map to the Pangenome with VG Giraffe [run_giraffe]
#### input

* sample [workflow inputs]
* fastqs [workflow inputs]
* gbz [workflow inputs]
* dist [workflow inputs]
* min [workflow inputs]
####output:
* gam

### VG Stats [vg_stats]
#### input
* gam [from run_giraffe]
####output:
* stats

### VG Surject GAM to BAM [run_surject]
#### input
* gam [from run_giraffe]
* sample [workflow inputs]
* library [workflow inputs + "-lib1"]
* gbz [workflow inputs]
#### output
* bam

### Samtools Sort BAM by Coordinates [samtools_sort]
The bam needs to be sorted by coordinate to call variants
#### input
* bam [bam from run_surject]
#### output
* sorted_bam

### Picard Markdup [markdup]
#### input
* bam [sorted_bam from samtools_sort]
#### output
* dedup_bam

### Samtools Stat [samtools_stat]
#### input
* bam [dedup_bam from markdup]
#### output
* stats [samtools stat file]

### Samtools Index
Deep variant and other tools require a BAI
#### input
* bam [dedup_bam from markdup]
#### output
* bai [bam index file]

### Deep Variant
* bam [dedup_bam from markdup]
* reference [workflow input]
#### output
* vcf

## Outputs
* gam [gam from run_giraffe]
* gam_stats [stats from run_stats]
* bam [dedup_bam from markdup]
* bai [bai from samtools_index]
* bam_stats [stats from samtools_stat]
* vcf [vcf from deep variant]
