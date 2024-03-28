# Minigraph-Cactus Graph Builder - MCGB

Build pangenome graphs with minigraph-cactus.

## Pipeline Chart
*Note: this workflow runs all 5 steps in one task*
```mermaid
  flowchart TB;
      i1([NAME]);
      i2([REF]);
      i3([SEQFILE])

      s1[CACTUS-MINIGRAPH];
      s2[CACTUS-GRAPHMAPH];
      s3[CACTUS-GRAPHMAP-SPLIT];
      s4[CACTUS-ALIGN];
      s5[CACTUS-GRAPHMAP-JOIN];

      o11([SV_GFA])
      o21([PAF])
      o22([FASTA])
      o31([CHROM_FILE])
      o32([CHROMS])
      o41([CHROMS_ALIGNMENTS])
      o42([VG])
      o51([GIRAFFE - DIST/GBZ/MIN])
      o52([HAL])

      o01([gfa])
      o03([stats])
      o04([vcf])
      o05([vcf_tbi])
      o06([raw_vcf])
      o07([raw_vcf_tbi])

      i1-->s1; i2-->s1; i3-->s1; s1-->o11;
      s1--WORDS-->s2;
      i1-->s2; s2-->o21;
      s2-->s3;
      i1-->s3; s3-->o31; s3-->o32
      s3-->s4;
      o31-->s4; o32-->s4; s4-->o41; s4-->o42;
      s4-->s5;
      o41-->s5; s5-->o51; s5-->o52
```

## Pipeline Files
* mcgb.wdl          - WDL pipeline
* mcgb.inputs.json  - pipeline inputs with place holders
* mcgb.outputs.yaml - steps and outputs to be copied after pipeline run
* mcgb.imports.zip  - imports used in the WDL
* mcgb.doc.md       - this file, documenting the pipeline

## Inputs
* name [String] - output base name
* ref [String] - the sequence name in the seqfile to as the "reference"
* seqfile [File] - tab sepqarated file of sequence names and URLs (files)
* docker [String] - docker to use
* cpu [Int] - cpus to request
* memory [Int] - memory to request

## Steps
### Map to the Pangenome with VG Giraffe [run_mcgb]
#### input
* name [workflow inputs]
* ref [workflow inputs]
* seqfile [workflow inputs]
#### output
* dist
* fasta
* gaf
* gbz
* gfa
* hal
* min
* paf
* stats
* sv_gfa
* vcf
* vcf_tbi
* raw_vcf
* raw_vcf_tbi
