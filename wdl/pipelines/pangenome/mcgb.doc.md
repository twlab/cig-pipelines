# Minigraph-Cactus Graph Builder - MCGB

Build pangenome graphs with minigraph-cactus.

## Pipeline Chart
```mermaid
  flowchart TB;
      i1([NAME]);
      i2([REF]);
      i3([SEQFILE])

      s1[[CACTUS-MINIGRAPH]];
      s2[[CACTUS-GRAPHMAPH]];
      s3[[CACTUS-GRAPHMAP-SPLIT]];
      s4[[CACTUS-ALIGN]];
      s5[[CACTUS-GRAPHMAP-JOIN]];

      o11([SV_GFA])
      o21([PAF])
      o22([FASTA])
      o23([GAF])
      o51([DIST/GBZ/MIN for GIRAFFE])
      o52([HAL])
      o53([VG])
      o54([GFA])
      o55([STATS])
      o56([VCF & TBI])
      o57([RAW VCF & TBI])

      i1-->s1; i2-->s1; i3-->s1; s1-->o11;
      o11-->s2; s2-->o21; s2-->o22;
      o21-->s3;
      s3--CHROMOSOMES-->s4;
      s4--CHROMOSOME ALIGNMENTs & HALs -->s5;
      s5-->o51; s5-->o52; s5-->o53; s5-->o54; s5-->o55; s5-->o56; s5-->o57
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
* docker [String] - docker to use (cactus)
* cpu [Int] - cpus to request
* memory [Int] - memory to request

## Steps
### run_cactus_minigraph
#### input
* name [workflow inputs]
* ref [workflow inputs]
* seqfile [workflow inputs]
#### output
* sv_gfa
### run_cactus_graphmap
#### input
* name [workflow inputs]
* ref [workflow inputs]
* seqfile [workflow inputs]
* sv_gfa [run_cactus_minigraph]
#### output
* fasta
* gaf
* paf
* seqfile [updated with minigraph fasta, not retained]
### run_cactus_graphmap_split
#### input
* ref [workflow inputs]
* seqfile [run_cactus_graphmap]
* paf [run_cactus_graphmap]
* sv_gfa [run_cactus_minigraph]
#### output
* chroms (not retained)
### run_cactus_align
#### input
* ref [workflow inputs]
* chroms [run_cactus_graphmap_split]
#### output
* alignments (not retained)
### run_cactus_graphmap_join
#### input
* name [workflow inputs]
* ref [workflow inputs]
* alignments [run_cactus_align]
#### output
* dist
* gbz
* gfa
* hal
* min
* stats
* vcf
* vcf_tbi
* raw_vcf
* raw_vcf_tbi
