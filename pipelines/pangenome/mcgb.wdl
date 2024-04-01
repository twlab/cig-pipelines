version development

import "wdl/structs/runenv.wdl"

workflow pangenome_mcgb {
  input {
    String name
    String ref
    File seqfile
    String docker = "mgibio/cactus:2.5.0-focal"
    Int cpu
    Int memory
  }

  RunEnv runenv_mcgb = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
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
    runenv=runenv_mcgb,
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
    Directory alignments
    RunEnv runenv
  }

  String jobstore = "jobstore"
  command <<<
    set -e
    export OMP_NUM_THREADS=1
    ln -s ~{alignments} alignments
    cactus-graphmap-join ~{jobstore} --vg alignments/*.vg --hal alignments/*.hal --outDir . --outName ~{name} --reference ~{ref} --vcf --giraffe clip --maxCores ~{runenv.cpu} --maxMemory ~{runenv.memory - 2}G --defaultDisk 100G --binariesMode local
  >>>

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
    #disks:  "local-disk ~{runenv.disks} SSD"
  }

  output {
    File dist = glob("*.dist")[0]
    File hal = glob("~{name}*.hal")[0]
    File gbz = glob("~{name}.gbz")[0]
    File gfa = glob("~{name}.gfa.gz")[0]
    File min = glob("~{name}.min")[0]
    File stats = glob("~{name}.stats.tgz")[0]
    File vcf = glob("~{name}.vcf.gz")[0]
    File vcf_tbi = glob("~{name}.vcf.gz.tbi")[0]
    File raw_vcf = glob("~{name}.raw.vcf.gz")[0]
    File raw_vcf_tbi = glob("~{name}.raw.vcf.gz.tbi")[0]
  }
}
