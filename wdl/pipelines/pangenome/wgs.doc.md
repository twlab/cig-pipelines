# Pangenome WGS Pipeline

Map WGS data to the pangenome graph using the "personalized pangenome method" using giraffe, polish alignemnt indels with freebayes left shift & abra2 realign indels, and call variants with haplotype caller and deep variant. Produces a GAM (w/ stats), positional sorted BAM (w/ BAI & stats), plus VCF (w/ TBI).

## Pipeline Chart
```mermaid
flowchart TB;
  i1([FASTQs]);
  i2([GBZ / HAP]);
  i3([SAMPLE])
  i4([REFERENCE NAME])

  s0[FastQC]
  s01[Generate KMERS]
  s1[VG GIRAFFE];
  s2[VG STATS];
  s3[EXTRACT REF FROM GRAPH];
  s4[VG SURJECT];
  s5[SAMTOOLS SORT];
  s6[FREEBAYES LEFT ALIGN];
  s7[SAMTOOLS INDEX];
  s8[GATK REALIGNER TARGET CREATOR];
  s9[BEDTOOLS EXPAND TARGETS];
  s10[ABRA2 REALIGN];
  s11[DEEPVARIANT];
  s12[SAMTOOLS STAT];

  o0([FastQC Output Files])
  o1([GAM])
  o2([GAM STAT])
  o3([BAM])
  o4([BAI])
  o5([BAM STAT])
  o6([VCF])
  o7([VCF TBI])

  i1---->s0;
  i1-->s01; i1-->s1; i2-->s1; i3-->s1; i3-->s4; i3-->s11;
  s01--KMERS-->s1; s1--GAM-->s2;
  s1--GAM-->s4--BAM-->s5--SORTED BAM-->s6--LEFT SHIFT BAM-->s8--TARGETS-->s9--EXPANDED TARGETS-->s10--BAM & BAI-->s11;
  s6--LEFT SHIFT BAM-->s10;
  s6--LEFT SHIFT BAM-->s7--BAI-->s10; s7--LEFT SHIFT BAI-->s8;
  s10--BAM-->s12;
  i2-->s3; i4-->s3; s3--EXTRACTED REF-->s10; s3--EXTRACTED REF-->s11;

  s0-->o0; s1-->o1; s2-->o2; s10-->o3; s10-->o4; s12-->o5; s11-->o6; s11-->o7;
```
## Pipeline Files
* wgs.wdl          - WDL pipeline
* wgs.inputs.json  - pipeline inputs with place holders
* wgs.outputs.yaml - steps and outputs to be copied after pipeline run
* wgs.imports.zip  - imports used in the WDL
* wgs.doc.md       - this file, documenting the pipeline

## Inputs
See [inputs json](https://github.com/twlab/cig-pipelines/blob/main/wdl/pipelines/pangenome/wgs.inputs.json) for all inputs.

* sample [String] - sample name for outputs
* fastqs [File] - an array of read1 and read2 fastqs 
* gbz [File] - giraffe pangenome GBZ
* hap [File] - vg haplotypes for the GBZ
* reference_name [String] - reference name to extract from the graph

## Outputs
See [inputs json](https://github.com/twlab/cig-pipelines/blob/main/wdl/pipelines/pangenome/wgs.outputs.yaml).
