# Pangenome Extract Reference FASTA Tool

Given a pangenome graph, extract the sequence of a reference.

## Pipeline Files
* extract-ref.wdl - WDL pipeline
* extract-ref.inputs.json - pipeline inputs with place holders
* extract-ref.outputs.yaml - steps and outputs to be copied after pipeline run
* extract-ref.imports.zip - imports used in the WDL
* extract-ref.doc.md - this file, documenting the pipeline

## Inputs
* name [String] - base name for outputs
* fastqs [File] - an array of an 2 arrays, one each for read1 and read2 fastqs
* idx [File] - tarred BWA index (made from build_idx workflow)
* known_sites - gzipped VCF of known sites like dbSNP

## Steps

### Run VG Paths and Samtools FAIDX & DICT
#### input
* name [inputs.name]
* gbz [inputs.gbz]
#### output
* fasta
* fai
* dict

## Outputs
* fasta [from run vg paths & samtools]
* fai [from run vg paths & samtools]
* dict [from run vg paths & samtools]
