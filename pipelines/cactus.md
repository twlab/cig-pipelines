# Cactus Pangeome Grpah Builder Pipeline

Given some sequence assemblies, build a pangenome grpah with cactus and minigraph.

## Pipeline Files
* cactus.wdl - WDL pipeline
* cactus.inputs.json - pipeline inputs with place holders
* cactus.outputs.yaml - steps and outputs to be copied after pipeline run
* cactus.imports.zip - imports used in the WDL
* cactus.imports.README - this file, documenting the pipeline

## Inputs
* sequences [File] - tab separated file of sequnce nams and locations
* reference_name [String] - the name of the sequence to use as the reference
* out_name [String] - base name for the outputs

### Resource Inputs
These resource inputs can also be given: docker, cpu, memory, and disks

## Steps
### Minigraph
Construct a minigraph in gfa format.
#### input
* sequences [workflow input]
* reference_name [workflow input]
* out_name [workflow input]
#### output

### Graphmap
Map each input assembly back to the graph with minigraph.
#### input
* sequences [workflow input]
* reference_name [workflow input]
* out_name [workflow input]
* gfa [output from minigraph]
#### output
* fasta [minigraph fasta]
* paf
* updated_sequences [sequences file to be copied and updated with local file locations]

### Graphmap Split
This step is optional, but helps in running by spliting the PAF into chromosomes.
#### input
* sequences [workflow input]
* reference_name [workflow input]
* out_name [workflow input]
* gfa [minigraph output]
* fasta [graphmap output]
* paf [graphmap output]
#### output
*  chroms [directory with seqfiles and pafs for each chromosome]

### Align
Compute the Cactus multiple genome alignment from the assembly-to-graph minigraph mappings.
#### input
* reference_name [workflow input]
* out_name [workflow input]
* fasta [graphmap output]
* chroms [graphmap split output]
#### output
* alignments [directory of chromosome alignments]

### Minigraph Join
Produce the final graph and indexes.
#### input
* alignments [align output]
* reference_name [workflow input]
* out_name [workflow input]
#### output
* dist [minigraphjoin output]
* gbz [minigraphjoin output]
* gfa [minigraphjoin output]
* hal [minigraphjoin output]
* min [minigraphjoin output]
* stats [minigraphjoin output]
* raw_vcf [minigraphjoin output]
* raw_vcf_tbi [minigraphjoin output]
* vcf [minigraphjoin output]
* vcf_tbi [minigraphjoin output]

## Outputs
* dist [minigraphjoin output]
* gbz [minigraphjoin output]
* gfa [minigraphjoin output]
* hal [minigraphjoin output]
* min [minigraphjoin output]
* stats [minigraphjoin output]
* raw_vcf [minigraphjoin output]
* raw_vcf_tbi [minigraphjoin output]
* vcf [minigraphjoin output]
* vcf_tbi [minigraphjoin output]
