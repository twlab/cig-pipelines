# Build Indexes for Bulk RNA Pipeline

Build STAR, RSEM, & KALLISTO indexes for the bulk RNA pipeline.

## Pipeline Chart
```mermaid
  flowchart TB;
      i1([REFERENCE]);
      i2([ANNOTATION]);
      i3([GENOME]);
      i4([ANNO_VERSION]);
      o1([TRANSCRIPTS FASTA]);
      o2([KALLISTO INDEX]);
      o3([STAR INDEX]);
      o4([RSEM INDEX]);
      s1[BUILD TRANSCRIPTS FASTA];
      s2[BUILD TRANSCRIPTS IDX];
      s3[BUILD STAR IDX];
      s4[BUILD RSEM IDX];
      i1-->s1;
      i2-->s1;
      s1--TRANSCRIPTS FASTA-->s2;
      i1-->s3;
      i2-->s3;
      i3-->s3;
      i4-->s3;
      i1-->s4;
      i2-->s4;
      i3-->s4;
      i4-->s4;
      s1-->o1;
      s2-->o2;
      s3-->o3;
      s4-->o4;
```

## Pipeline Files
* build-idx.wdl - WDL pipeline
* build-idx.inputs.json - pipeline inputs with place holders
* build-idx.outputs.yaml - steps and outputs to be copied after pipeline run
* build-idx.imports.zip - imports used in the WDL
* build-idx.doc.md - this file, documenting the pipeline

## Inputs
* reference [File] - sequences FASTA GZ
* annotation [File] - annotations GTF GZ
* genome [String] - genome name (GRCh38)
* anno_version - annotation version (gencode_v36)

## Steps
### Build Transcripts FASTA [build_transripts_fasta]
#### input
* reference [inputs.reference]
* annotation [inputs.annotation]
#### output
* transcripts [File] - transcript sequences FASTA

## Build Transcripts (Kallisto) Index [build_transcripts_index]
####input
* reference [build_transcripts_fasta.transcripts]
####output
* index [File] - kallisto index

### Build STAR Index [build_star_index]
#### input
* reference [inputs.reference]
* annotation [inputs.annotation]
* genome [inputs.genome]
* anno_version [inputs.anno_version]
#### output
* index [File] - TAR with STAR index plus support files

### Build RSEM Index [build_rsem_index]
#### input
* reference [inputs.reference]
* annotation [inputs.annotation]
* genome [inputs.genome]
* anno_version [inputs.anno_version]
#### output
* index - RSEM index

## Outputs
* kallisto_index [build_transcripts_index.index]
* rsem_index [build_rsem_index.index]
* star_index [build_star_index.index]
* trancripts [build_transcripts.transcripts]
