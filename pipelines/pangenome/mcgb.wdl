version development

import "wdl/structs/runenv.wdl"

workflow pangenome_mcgb {

  input {
    String name
    String ref
    File seqfile
    String docker = "ebelter/cactus:2.5.0-20.04"
    Int cpu
    Int memory
  }

  RunEnv runenv_mcgb = {
    "docker": docker,
    "cpu": cpu,
    "memory": memory,
    "disks": 20,
  }

  call run_mcgb { input:
    name=name,
    ref=ref,
    seqfile=seqfile,
    runenv=runenv_mcgb,
  }

  output {
    File gfa = run_mcgb.gfa
    File paf = run_mcgb.paf
    File fasta = run_mcgb.fasta
  }
}

task run_mcgb {
  input {
    String name
    String ref
    File seqfile
    RunEnv runenv
  }

  String jobstore = "jobstore"
  String sv_gfa = "~{name}.sv.gfa" # GFA="${OUTDIR}/${NAME}.sv.gfa"
  String paf = "~{name}.paf" # PAF="${OUTDIR}/${NAME}.paf"
  String fasta = "~{name}.sv.gfa.fasta" # FASTA="${OUTDIR}/${NAME}.sv.gfa.fa"
  String chroms = "chroms"
  String alignments = "chrom-alignments"
  command <<<
    set -e
    cactus-minigraph ~{jobstore} ~{seqfile} ~{sv_gfa} --reference ~{ref} --maxCores ~{runenv.cpu} --maxMemory ~{runenv.memory - 2}G --defaultDisk 100G --binariesMode local
    
    cactus-graphmap ~{jobstore} ~{seqfile} ~{sv_gfa} ~{paf} --outputFasta ~{fasta} --reference ~{ref} --maxMemory ~{runenv.memory - 2}G --defaultDisk 100G --binariesMode local
    
    cactus-graphmap-split ~{jobstore} ~{seqfile} ~{sv_gfa} ~{paf} --outDir ~{chroms} --reference ~{ref} --maxMemory ~{runenv.memory - 2}G --defaultDisk 100G --binariesMode local
    
    chromfile=$(find ~{chroms} -name chromfile.txt)
    cactus-align ~{jobstore} $chromfile ~{alignments} --batch --pangenome --reference ~{ref} --outVG --maxMemory ~{runenv.memory - 2}G --defaultDisk 100G --binariesMode local
    
    cactus-graphmap-join ~{jobstore} --vg ~{alignments}/*.vg --hal ~{alignments}/*.hal --outDir . --outName ~{name} --reference ~{ref} --vcf --giraffe clip --maxMemory ~{runenv.memory - 2}G --defaultDisk 100G --binariesMode local
  >>>

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
    #disks:  "local-disk ~{runenv.disks} SSD"
  }

  output {
    File dist = glob("*.dist")[0]
    File fasta = glob("~{fasta}")[0]
    File gaf = glob("~{name}.gaf.gz")[0]
    File gbz = glob("~{name}.gbz")[0]
    File gfa = glob("~{name}.gfa.gz")[0]
    File hal = glob("~{name}*.hal")[0]
    File min = glob("~{name}.min")[0]
    File paf = glob("~{paf}")[0]
    File stats = glob("~{name}.stats.tgz")[0]
    File sv_gfa = glob("~{sv_gfa}")[0]
    File vcf = glob("~{name}.vcf.gz")[0]
    File vcf_tbi = glob("~{name}.vcf.gz.tbi")[0]
    File raw_vcf = glob("~{name}.raw.vcf.gz")[0]
    File raw_vcf_tbi = glob("~{name}.raw.vcf.gz.tbi")[0]
  }
}
