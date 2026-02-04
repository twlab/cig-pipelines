# Distortion Map Pipeline


## Pipeline Chart
```mermaid
---
title: Distortion Map Pipeline
---
flowchart TB;
  i1([Query Chromsomes TSV]);

  s21[[Generate Simulated Reads]];
  s22[[Extract SOURCE Read Poistions]];
  s23[[Liftover SOURCE Read Poisitons to REF]];

  s31[[Map Reads to REF]];
  s32[[BAM to BED]];

  s41[[Map Reads to Query]];
  s42[[BAM to BED]];
  s43[[Liftover Alignments to REF]];

  s61[[Integrate Data into SQL Database]];
  s62[[Create Intervals]]
  s63[[Generate Count Matrices]];

  s101[[Merge & Normalize REF Count Matrices]];
  s102[[Merge & Normalize SOURCE Count Matrices]];
  s111[[Merge Intervals **NEW**]];
  s121[[Calculate Distortion Metrics]];

  o1([Metrics]);
 
  i1-->s21;
  subgraph sg1 ["SCATTER BY CHROMOSOME"]
    s21--Simulated Reads-->s31;
    s21--Simulated Reads-->s41;
    s21--Simulated Reads-->s22-->s23--Simulated Reads Lift Over Positions-->s61;
    s22--Simulated Reads Positions-->s61
    subgraph sgr ["QUERY"]
      s41--Alignments-->s42--BED-->s43
      s42--Query Alignment Poisitions-->s61;
    end
    subgraph sgq ["REFERENCE"]
      s31--Alignments-->s32
    end
    s32--Ref Alignment Positions-->s61;
    s43--Query Liftover Positions-->s61;
    s61--DB-->s62--DB / Intervals-->s63;
  end
  s63--REF Matrices per Chromosome-->s101--Merged REF Matrices-->s121
  s63--Intervals per Chromosome-->s111--Merged Intervals-->s121
  s63--SOURCE Matrices per Chromosome-->s102--Merged SOURCE Matrices-->s121
  s121-->o1
```
