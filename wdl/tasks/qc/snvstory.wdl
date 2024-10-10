version development

import "../structs/runenv.wdl"

task run_igm_churchill_ancestry {
  input {
    File input_vcf
    Directory resource
    String genome_ver = "38"
    String mode = "WGS"
    RunEnv runenv
  }

  String output_dir = "out"
  command <<<
    set -x
    mkdir -p ~{output_dir}
    export TMP_DIR=$(readlink -f ~{output_dir})
    export NUMBA_CACHE_DIR="${TMP_DIR}/numba"
    mkdir -p "${NUMBA_CACHE_DIR}"
    /opt/conda/bin/python3 -m igm_churchill_ancestry \
      --path ~{input_vcf} \
      --resource ~{resource} \
      --output-dir ~{output_dir} \
      --genome-ver ~{genome_ver} \
      --mode ~{mode}
  >>>

  output { # *_1kGP_umap.html *_genotyping.vcf_PAN001.csv *_gnomAD_umap.html *.pdf *_SGDP_umap.html
    File csv_file = glob("~{output_dir}/output/*.csv")[0]
    File pickle_file = glob("~{output_dir}/output/*.pickle")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "${runenv.memory} GB"
    #disks: runenv.disks
  }
}
