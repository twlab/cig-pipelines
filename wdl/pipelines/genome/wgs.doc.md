# Pangenome WGS Pipeline

Map WGS data to the linear geome with BWA, optionally left shift/realign bam, and call variants with deep variant. Produces sorted BAM w/ BAI, VCF, and stats.

## Pipeline Chart
```mermaid
  flowchart TB;
      i1([SAMPLE])
      i2([FASTQs]);
      i3([BWA IDX]);
      i4([TARGETS_EXPANSION_BASES])

      s1[UNTAR IDX];
      s2[BWA ALIGN];
      s3[SAMTOOLS SORT];
      s4[FREEBAYES LEFT ALIGN];
      s5[SAMTOOLS INDEX];
      s6[GATK REALIGNER TARGET CREATOR];
      s7[BEDTOOLS EXPAND TARGETS];
      s8[ABRA2 REALIGN];
      s9[DEEPVARIANT];
      s10[SAMTOOLS STAT];

      o1([BAM])
      o2([BAI])
      o3([BAM STAT])
      o4([VCF])
      o5([VCF TBI])

      i3-->s1;
      s1--REF PATH-->s2; i1-->s2; i2-->s2;
      s2--BAM-->s3;
      s3--SORTED BAM-->s4;
      s4--LEFT SHIFT BAM-->s5;
      s4--LEFT SHIFT BAM-->s6; s5--LEFT SHIFT BAI-->s6;
      s6--TARGETS-->s7; i4-->s7;
      s4--LEFT SHIFT BAM-->s8; s7--EXPANDED TARGETS-->s8;
      s8--REALIGN BAM/BAI-->s9;
      s8--REALIGN BAM-->s10;
      
      s8-->o1; s8-->o2; s9-->o4; s9-->o5; s10-->o3;
```
## Pipeline Files
* wgs.wdl          - WDL pipeline
* wgs.inputs.json  - pipeline inputs with place holders
* wgs.outputs.yaml - steps and outputs to be copied after pipeline run
* wgs.doc.md       - this file, documenting the pipeline

## Inputs
* sample [String] - sample name for outputs
* fastqs [File] - an array of read1 and read2 fastqs 
* idx [File] - bwa index with reference fastqa, fai, and dict

