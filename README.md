# CIG@MGI Pipelines

The Collaborative and Intergrative Genomics (CIG) is a group in the [McDonnell Genome Institute](https://www.genome.wustl.edu/) (MGI) at the [Washington University School of Medicine](https://medicine.wustl.edu/) (WUSM).

## Overview

In this repo, we share genomics pipelines, tools, and tasks, mostly in the form of [WDL](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md). Corresponding Dockerfiles and scripts are at [CIG Containters](https://github.com/twlab/cig-containers).

## Repo Structure

| Path          | Description |
| ---           | --- |
| scripts       | helpers scripts for pipelines; these are not used in pipelines/tools |
| wdl/pipelines | full start to end workflows that produce outputs from many tasks |
| wdl/tasks     | wrapped command line interfaces and scripts (must be incorporated into tools/pipelines) |
| wdl/tools     | stand alone workflows that combine a few tasks to produce outputs |
| wdl/structs   | data structures and types for pipelines, tools, and tasks |

# Building WDL Imports
WDLs have imports that bring in functionality from other WDL files. When running **cromwell**, a zip file of the imports is required. These WDL imports ZIP files are not included in the repo, but can be easily built with the *build_imports* script. From the main directory run the *build_imports* script giving the WDL file as the sole input. The output will be next to the WDL, named with the basename of the WDL file with the extension *.wdl* replaced with *.imports.zip*. Here is an example building the imports for *wdl/tools/pacbio/align_output_cram.wdl*:

```
$ cd cig-pipelines
$ bash ./scripts/build_imports wdl/tools/pacbio/align_output_cram.wdl 
WDL: wdl/tools/pacbio/align_output_cram.wdl
Archive:  align_output_cram.imports.zip
  Length      Date    Time    Name
---------  ---------- -----   ----
        0  2025-08-13 14:14   wdl/
        0  2025-08-13 14:14   wdl/structs/
      285  2025-08-13 14:14   wdl/structs/runenv.wdl
        0  2025-08-13 14:14   wdl/tasks/
        0  2025-08-13 14:14   wdl/tasks/pacbio/
     1568  2025-08-13 14:14   wdl/tasks/pacbio/pbmm2.wdl
---------                     -------
     1853                     6 files
Created imports file: wdl/tools/pacbio/align_output_cram.imports.zip

```

The *tools* and *pipelines* imports include the absolute paths of import WDL files. The *tasks* imports, which are typically limited, include the relative path of the imports. Pipeline and tool WDLs can use task WDLs, and task WDLs should only use other task WDLs.

