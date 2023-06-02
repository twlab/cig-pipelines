# Align and Call Pipeline

Align WGS and call variants.  This pipeline is a standard for aligning WGS data to a reference using bwa-mem, then calling variants with GATK haplocaller.

## Files
* align-and-call.wdl - WDL pipeline
* align-and-call.inputs.json - pipeline inputs with place holders
* align-and-call.outputs.yaml - steps and outputs to be copied after pipeline run
* align-and-call.imports.zip - imports used in the WDL

## Inputs

* name [String] - base name for outputs
* fastqs [File] - read fastqs
* idx [File] - tarred BWA index (made from build_idx workflow)
* known_sites - gzipped VCF of known sites like dbSNP

## Steps

### Untar the BWA Reference
#### input
* idx [inputs.idx]
#### output
* reference - untarred BWA idx

### Align with BWA MEM
#### input
* name [workflow inputs]
* fastqs [workflow inputs]
* reference [output from untar idx]
####output:
* bam

### Collect Samtools Stat
#### input
* bam [output from align]
#### output
* stats [samtools stats file]

### Sort the BAM by Coordinates
The bam needs to be sorted by coordinate to call variants
#### input
* bam [output from align]
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
