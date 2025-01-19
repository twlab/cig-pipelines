# Pangenome WGS Pipeline

Map WGS data to the linear geome with BWA, optionally left shift/realign bam, and call variants with deep variant. Produces sorted BAM w/ BAI, VCF, and stats.

## Pipeline Chart
```mermaid
flowchart TB;
  i1([SAMPLE])
  i2([FASTQs]);
  i3([BWA IDX]);
  i4([TARGETS_EXPANSION_BASES])
  s0[[FastQC]];
  s1[[UNTAR IDX]];
  s2[[BWA ALIGN]];
  s3[[SAMTOOLS SORT]];
  s4[[FREEBAYES LEFT ALIGN]];
  s5[[SAMTOOLS INDEX]];
  s6[[GATK REALIGNER TARGET CREATOR]];
  s7[[BEDTOOLS EXPAND TARGETS]];
  s8[[ABRA2 REALIGN]];
  s9[[PICARD MARKDUP]];
  s10[[SAMTOOLS INDEX]];
  s11[[DEEPVARIANT]];
  s12[[SAMTOOLS STAT]];

  i2-->s0;
  i3-->s1;
  s1--REF PATH-->s2; i1-->s2; i2-->s2;
  s2--BAM-->s3;
  s3--SORTED BAM-->s4;

  subgraph sg1 [Left Align/Realign Bam]
    s4--LEFT SHIFT BAM-->s5;
    s4--LEFT SHIFT BAM-->s6; s5--LEFT SHIFT BAI-->s6;
    s6--TARGETS-->s7; i4-->s7;
    s4--LEFT SHIFT BAM-->s8; s7--EXPANDED TARGETS-->s8;
  end;
  s8--REALIGN BAM-->s9;
  s9--REALIGN/MARKDUP BAM-->s10;
  s9--REALIGN/MARKDUP BAM-->s11;
  s9--REALIGN/MARKDUP BAM-->s12;
  s10--BAI-->s11;

  o0([FastQC Output Files])
  o1([BAM])
  o2([BAI])
  o3([BAM STAT])
  o4([VCF])
  o5([VCF TBI])

  s0-->o0; s9-->o1; s10-->o2; s11-->o4; s1-->o5; s12-->o3;
```
## Pipeline Files
* wgs.wdl          - WDL pipeline
* wgs.inputs.json  - pipeline inputs with place holders
* wgs.outputs.yaml - steps and outputs to be copied after pipeline run
* wgs.doc.md       - this file, documenting the pipeline

## Inputs
See [inputs json](https://github.com/twlab/cig-pipelines/blob/main/wdl/pipelines/genome/wgs.inputs.json) for all inputs.

* sample [String] - sample name for outputs
* fastqs [File] - an array of read1 and read2 fastqs 
* idx [File] - tarred file with bwa index, reference fastqa, fai, and dict

## Outputs
See [inputs json](https://github.com/twlab/cig-pipelines/blob/main/wdl/pipelines/genome/wgs.outputs.yaml).
 
