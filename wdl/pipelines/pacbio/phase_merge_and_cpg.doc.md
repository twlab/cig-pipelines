# SMaHT Long Read QC Pipeline

PacBio long read QC by aligning nad calling variants with DeepVariant. Additionally run CPG and  Qualimap. Use these VCFs to determine sex and ancestory. These results are collated into a table report.

## Pipeline Chart
```mermaid
flowchart LR
  i1([LIBINFO TSV]);
  i2([REF MMI]);
  i3([REF FASTA]);
  i4([REF FAI]);
  i51([SNVSTORY RESOURCE]);
  i61([CENTRIFUGER DBS]);
  i63([VERIFYBAMID RESOURCE]);

  s1[[PBMM2 ALIGN]];
  s2[[DEEP VARIANT]];
  s31[[QUALIMAP]];
  s32[[CPG]];
  s41[[SAMTOOLS STATS]];
  s51[[SNVSTORY ANCESTRY]];
  s52[[CHECK SEX]];
  s61[[DOWNSAMPLE & CENTRIFGUER]];
  s62[[HAPLOCHECK]];
  s63[[VERIFYBAMID]];
  s71[GENERATE QC REPORT];

  o11{{BAMs & BAIs}};
  o21{{VCFs & TBIs plus HISTOs per Library/FASTQs}};
  o31{{RESULTs & REPORTs}};
  o32{{BEDs & BIGWIGs}};
  o71{{QC REPORT}};

  i61-->s61;
  i63-->s63;
  subgraph sg1 [SCATTER BY Library/BAMs]
    s1--BAM/BAI-->s2;
    s1--BAM-->s31;
    s1--BAM/BAI-->s32;
    s1--BAM-->s41;
    s1--BAM-->s61;
    s1--BAM/BAI-->s62;
    s1--BAM/BAI-->s63;
    s2--VCF-->s51;
    s2--VCF-->s52;

  end
  s41--STATS-->s71;
  s51--TSVs-->s71;
  s52--DATAFRAMEs-->s71;
  s61--CONTIMATION REPORT-->s71;
  s62--CONTIMATION REPORT-->s71;
  s63--CONTIMATION REPORT-->s71;
  i1-->s1; i2-->s1;
  i3-->s2; i4-->s2;
  i51-->s51;
  s1-->o11;
  s2-->o21;
  s31-->o31;
  s32-->o32;
  s71-->o71;
```

## Pipeline Files
* long_read_qc.wdl - WDL pipeline
* long_read_qc.inputs.json - pipeline inputs with place holders
* long_read_qc.outputs.yaml - steps and outputs to be copied after pipeline run
* long_read_qc.doc.md - this file, documenting the pipeline

## Imports
To generate imports, use /app/scritps/build_imports in smaht-qc docker image. See main repo README for example.

## Inputs
See *long_read_qc_batch.inputs.json* file. To generate the libinfo TSV with the smaht command. See main repo README for example.

## Outputs
See *long_read_qc_batch.outputs.yaml* file.
