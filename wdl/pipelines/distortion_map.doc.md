# Distortion Map Pipeline


## Pipeline Chart
```mermaid
---
title: Distortion Map Pipeline
---
flowchart TB;
  s11[[SCATTER BY CHOMOSOME Extract Chromosome FASTAs]]
  s12[[MiniMap2 QUERY to REF]]

  s21[[Generate Simulated Reads]];
  s22[[Extract SOURCE Read Poistions]];
  s23[[Liftover SOURCE Read Poisitons to REF]];

  s30[[Build REF BWA Chromsome Index]]
  s31[[Map Reads to REF]];
  s32[[BAM to BED]];

  s40[[Build QUERY BWA Chromsome Index]]
  s41[[Map Reads to Query]];
  s42[[BAM to BED]];
  s43[[Liftover Alignments to REF]];

  s61[[Integrate Data into SQL Database]];
  s62[[Create Intervals]]
  s63[[Generate Count Matrices]];

  s101[[Merge & Normalize REF Count Matrices]];
  s102[[Merge & Normalize SOURCE Count Matrices]];
  s111[[Merge Intervals]];
  s121[[Calculate Distortion Metrics]];

  o1([Metrics]);
 
  subgraph sg1 [" "]
    subgraph qry ["QUERY"]
      s40--QUERY BWA IDX-->s41--Alignments-->s42--BED-->s43;
    end

    subgraph sim ["SIMULATE READS"]
      s21--Simulated Reads-->s22
      s22--Simulated Reads-->s23;
    end

    s11--REF & QUERY CHR-->s12;
    s11--QUERY CHR-->s21;
    s11--REF CHR-->s30;
    s11--QUERY CHR-->s40;
    s12--PAF-->s23;
    s12--PAF-->s43;

    s21--Simulated Reads-->s31;
    s21--Simulated Reads-->s41;
    s22--Simulated Reads Positions-->s61

    subgraph ref ["REFERENCE"]
      s30--REF BWA IDX-->s31--Alignments-->s32;
    end

    s23--Simulated Reads Lift Over Positions-->s61
    s32--Ref Alignment Positions-->s61;
    s43--Query Liftover Positions-->s61;
    s61--DB-->s62--DB / Intervals-->s63;
  end
  s63--REF Matrices per Chromosome-->s101--Merged REF Matrices-->s121;
  s63--Intervals per Chromosome-->s111--Merged Intervals-->s121;
  s63--SOURCE Matrices per Chromosome-->s102--Merged SOURCE Matrices-->s121;
  s121-->o1
```
