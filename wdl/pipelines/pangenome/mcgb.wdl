version development

import "wdl/structs/runenv.wdl"
import "wdl/tasks/pangenome/panacus.wdl"

workflow pangenome_mcgb {
  input {
    String name
    String ref
    File seqfile
    Array[String] graph_types = ["clip", "full"]
    String docker = "mgibio/cactus:2.5.0-focal"
    Int cpu
    Int memory
  }

  # Run Envss in order of use
  RunEnv runenv_mcgb = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
    "disks": 20,
  }

  RunEnv runenv_panacus = {
    "docker": "mgibio/panacus:0.2.3-buster",
    "cpu": 1,
    "memory": 24,
    "disks": 20,
  }

  call run_cactus_minigraph { input:
    name=name,
    ref=ref,
    seqfile=seqfile,
    runenv=runenv_mcgb,
  }

  call run_cactus_graphmap { input:
    name=name,
    ref=ref,
    seqfile=seqfile,
    sv_gfa=run_cactus_minigraph.sv_gfa,
    runenv=runenv_mcgb,
  }

  call run_cactus_graphmap_split { input:
    ref=ref,
    seqfile=run_cactus_graphmap.updated_seqfile,
    minigraph_fasta=run_cactus_graphmap.minigraph_fasta,
    paf=run_cactus_graphmap.paf,
    sv_gfa=run_cactus_minigraph.sv_gfa,
    runenv=runenv_mcgb,
  }

  call run_cactus_align { input:
    ref=ref,
    chroms=run_cactus_graphmap_split.chroms,
    runenv=runenv_mcgb,
  }

  call run_cactus_graphmap_join { input:
    name=name,
    ref=ref,
    alignments=run_cactus_align.alignments,
    graph_types=graph_types,
    runenv=runenv_mcgb,
  }

  scatter (gfa_gz in run_cactus_graphmap_join.gfa) {
    call panacus.run_panacus_hist as panacus { input:
      gfa_gz=gfa_gz,
      runenv=runenv_panacus,
    }
  }
}

task run_cactus_minigraph {
  input {
    String name
    String ref
    File seqfile
    RunEnv runenv
  }

  String jobstore = "jobstore"
  String sv_gfa = "~{name}.sv.gfa" # GFA="${OUTDIR}/${NAME}.sv.gfa"
  command <<<
    set -e
    export OMP_NUM_THREADS=1
    cactus-minigraph ~{jobstore} ~{seqfile} ~{sv_gfa} --reference ~{ref} --maxCores ~{runenv.cpu} --maxMemory ~{runenv.memory - 2}G --defaultDisk 100G --binariesMode local
  >>>

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
    #disks:  "local-disk ~{runenv.disks} SSD"
  }

  output {
    File sv_gfa = glob("~{sv_gfa}")[0]
  }
}

task run_cactus_graphmap {
  input {
    String name
    String ref
    File seqfile
    File sv_gfa
    RunEnv runenv
  }

  String jobstore = "jobstore"
  String paf = "~{name}.paf" # PAF="${OUTDIR}/${NAME}.paf"
  String fasta = "~{name}.sv.gfa.fasta" # FASTA="${OUTDIR}/${NAME}.sv.gfa.fa"
  String updated_seqfile = basename(seqfile)
  command <<<
    set -e
    export OMP_NUM_THREADS=1
    cp ~{seqfile} ~{updated_seqfile}
    cactus-graphmap ~{jobstore} ~{updated_seqfile} ~{sv_gfa} ~{paf} --outputFasta ~{fasta} --reference ~{ref} --maxCores ~{runenv.cpu} --maxMemory ~{runenv.memory - 1}G --defaultDisk 100G --binariesMode local
  >>>

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
    #disks:  "local-disk ~{runenv.disks} SSD"
  }

  output {
    File minigraph_fasta = glob("~{fasta}")[0]
    File gaf = glob("*.gaf.gz")[0]
    File paf = glob("~{paf}")[0]
    File updated_seqfile = glob("~{updated_seqfile}")[0]
  }
}

task run_cactus_graphmap_split {
  input {
    String ref
    File seqfile
    File sv_gfa
    File paf
    File minigraph_fasta
    RunEnv runenv
  }

  String jobstore = "jobstore"
  command <<<
    set -e
    export OMP_NUM_THREADS=1
    sed -i 's,_MINIGRAPH_\t.*,_MINIGRAPH_\t~{minigraph_fasta},' ~{seqfile}
    cactus-graphmap-split ~{jobstore} ~{seqfile} ~{sv_gfa} ~{paf} --outDir chroms --reference ~{ref} --maxCores ~{runenv.cpu} --maxMemory ~{runenv.memory - 1}G --defaultDisk 100G --binariesMode local
  >>>

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
    #disks:  "local-disk ~{runenv.disks} SSD"
  }

  output {
    Directory chroms = "chroms"
  }
}

task run_cactus_align {
  input {
    String ref
    Directory chroms
    RunEnv runenv
  }

  String jobstore = "jobstore"
  command <<<
    set -e
    export OMP_NUM_THREADS=1
    ln -s ~{chroms} chroms
    find chroms/seqfiles/ -type f | while read -r sf; do sed -i 's,file:///.*/execution/,,' $sf; done
    cactus-align ~{jobstore} chroms/chromfile.txt alignments --batch --pangenome --reference ~{ref} --outVG --maxCores ~{runenv.cpu} --maxMemory ~{runenv.memory - 1}G --defaultDisk 100G --binariesMode local
  >>>

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
    #disks:  "local-disk ~{runenv.disks} SSD"
  }

  output {
    Directory alignments = "alignments"
  }
}

task run_cactus_graphmap_join {
  input {
    String name
    String ref
    Array[String] graph_types
    Directory alignments
    RunEnv runenv
  }

  String jobstore = "jobstore"
  # --giraffe [GIRAFFE [GIRAFFE ...]]
  #           Generate Giraffe (.dist, .min) indexes for the given graph type(s). Valid types are 'full', 'clip' and 'filter'. If
  #           not type specified, 'filter' will be used (will fall back to 'clip' than full if filtering, clipping disabled,
  #           respectively). Multiple types can be provided seperated by a space
  command <<<
    set -e
    export OMP_NUM_THREADS=1
    ln -s ~{alignments} alignments
    cactus-graphmap-join ~{jobstore} --vg alignments/*.vg --hal alignments/*.hal --outDir . --outName ~{name} --reference ~{ref} --vcf ~{sep=' ' graph_types} --gfa ~{sep=' ' graph_types} --giraffe ~{sep=' ' graph_types} --maxCores ~{runenv.cpu} --maxMemory ~{runenv.memory - 2}G --defaultDisk 100G --binariesMode local
  >>>

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
    #disks:  "local-disk ~{runenv.disks} SSD"
  }

  output {
    Array[File] dist = glob("*.dist")
    Array[File] hal = glob("*.hal")
    Array[File] gbz = glob("*.gbz")
    Array[File] gfa = glob("*.gfa.gz")
    Array[File] min = glob("*.min")
    Array[File] vcf = glob("*.vcf.gz")
    Array[File] vcf_tbi = glob("*.vcf.gz.tbi")
    File stats = glob("*.stats.tgz")[0]
  }
}
