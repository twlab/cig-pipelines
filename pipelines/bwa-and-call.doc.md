# BWA and Call Pipeline

Map WGS fastqs using BWA-MEM and call variants with deep variant. This is the MGI Wang Lab standard for processing WGS data.

## Pipeline Files
* bwa-and-call.wdl - WDL pipeline
* bwa-and-call.inputs.json - pipeline inputs with place holders
* bwa-and-call.outputs.yaml - steps and outputs to be copied after pipeline run
* bwa-and-call.imports.zip - imports used in the WDL
* bwa-and-call.imports.README - this file, documenting the pipeline

## Inputs
* name [String] - base name for outputs
* fastqs [File] - an array of an 2 arrays, one each for read1 and read2 fastqs
* idx [File] - tarred BWA index (made from build_idx workflow)

## Steps
### Untar the BWA Reference
#### input
* idx [inputs.idx]
#### output
* reference - untarred BWA idx

### bwa with BWA MEM
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

### Samtools Sort BAM by Coordinates
The bam needs to be sorted by coordinate to call variants
#### input
* bam [output from samtools merge]
#### output
* sorted_bam

### Samtools Stat
#### input
* bam [output from align]
#### output
* stats [samtools stat file]

### Picard Mark Duplicates
#### input
* name [workflow inputs]
* bam [output bam from samtools sort]
#### output
* dedup_bam
* metrics

### Picard BQSR
[GATK BQSR doc](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-)
#### input
* bam [output bam from picard markdup]
* reference [output from untar idx]
* known_sites [workflow input]
#### output
* recalibrated bam
 
### Samtools Index
#### input
* bam [output from samtools merge]
#### output
* bai [bam index file]

### Deep Variant
* bam [otuput bam from gatk bqsr]
* reference [output from untar idx]
#### output
* vcf

## Outputs
* bam [from gatk_bqsr.recal_bam]
* bai [form samtools index]
* vcf [from deep variant]
* stats [from samtools stat]
* dedup_metrics [from markdup]
