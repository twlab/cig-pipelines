# CIG@MGI Pipelines

The Collaborative and Intergrative Genomics (CIG) is a group in the [McDonnell Genome Institute](https://www.genome.wustl.edu/) (MGI) at the [Washington University School of Medicine](https://medicine.wustl.edu/) (WUSM).

## Overview

In this repo, we share genomics pipelines, tools, and tasks, mostly in the form of [WDL](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md). Corresponding Dockerfiles and scripts are at [CIG Containters](https://github.com/twlab/cig-containers).

## Repo Structure

| Path      | Description |
| ---       | --- |
| pipelines | full start to end workflows that produce outputs from many tasks |
| tasks     | wrapped command line interfaces and scripts (must be incorporated into tools/pipelines) |
| tools     | stand alone workflows that combine a few tasks to produce outputs |
| structs   | data structures and types for pipelines, tools, and tasks |
| import    | assorted unassimilated pipelines, code, etc. |
