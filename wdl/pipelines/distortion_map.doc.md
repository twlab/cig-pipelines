# Distortion Map Pipeline


## Pipeline Chart
```mermaid
---
title: Distortion Map Pipeline
---
flowchart TB;
  i1([Query FASTA]);

  s1[[Split Query by Chromosome]];
  s2[[Generate Simulated Reads]];
  s31[[Map Reads to Query]];
  s32[[Map Reads to Ref]];
  s41[[BAM to BED]];
  s42[[BAM to BED]];
  s51[[Liftover]];
  s52[[Liftover]];
  s61[[Integrate Data into SQL Database]];
  s62[[Integrate Data into SQL Database]];
  s71[[Generate & Process Intervals]];
  s72[[Generate & Process Intervals]];
  s8[[Generate & Normalize Count Matrices]];
  s9[[Calculate Distortion Metrics]];

  o1([Graph])

  i1-->s1;
  s1--Query Chromosome FASTAs-->s2;
  subgraph sg1 ["SCATTER BY CHROMOSOME"]
    s2-->s31;
    s2-->s32;
    subgraph sgq ["QUERY"]
      s31-->s41-->s51-->s61-->s71;
    end
    subgraph sgr ["REF"]
      s32-->s42-->s52-->s62-->s72;
    end
  end
  s71-->s8; s72-->s8;
  s8-->s9;

  s9-->o1
```
