# Pangenome WGS Pipeline

Map WGS data to the pangenome graph using giraffe, and call variants with haplotype caller and deep variant. Produces a GAM, sorted BAM & BAI, plus VCF & stats.

## Pipeline Chart
```mermaid
  flowchart TB;
      i1([FASTQs]);
      i2([GBZ / DIST / MIN]);
      i3([SAMPLE])
      i4([REFERENCE NAME])

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

      o1([GAM])
      o2([GAM STAT])
      o3([BAM])
      o4([BAI])
      o5([BAM STAT])
      o6([VCF])
      o7([VCF TBI])

      i1-->s1; i2-->s1; i3-->s1; i3-->s4; i3-->s11;
      s1--GAM-->s2;
      s1--GAM-->s4--BAM-->s5--SORTED BAM-->s6--LEFT SHIFT BAM-->s8--TARGETS-->s9--EXPANDED TARGETS-->s10--BAM & BAI-->s11;
      s6--LEFT SHIFT BAM-->s10;
      s6--BAM-->s7--BAI-->s10; s7--LEFT SHIFT BAI-->s8;
      s10--BAM-->s12;
      i2-->s3; i4-->s3; s3--EXTRACTED REF-->s10; s3--EXTRACTED REF-->s11;
      s1-->o1; s2-->o2; s10-->o3; s10-->o4; s12-->o5; s11-->o6; s11-->o7;
```
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
* dist [File] - pangenome dist
* min [File] - pangenome min
* sample [String] - reference name to extract from the graph

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
