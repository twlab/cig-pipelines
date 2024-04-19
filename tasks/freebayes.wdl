version development

import "../structs/runenv.wdl"

task run_left_shift_bam {
  input {
    File in_bam_file
    File in_reference_file
    File in_reference_index_file
    RunEnv runenv
  }

  String out_prefix = basename(in_bam_file, ".bam")
  command <<<
    # Set the exit code of a pipeline to that of the rightmost command 
    # to exit with a non-zero status, or zero if all commands of the pipeline exit 
    set -o pipefail
    # cause a bash script to exit immediately when a command fails 
    set -e
    # cause the bash shell to treat unset variables as an error and exit immediately 
    set -u
    # echo each line of the script to stdout so we can see what is happening 
    set -o xtrace
    #to turn off echo do 'set +o xtrace' 

    # Reference and its index must be adjacent and not at arbitrary paths
    # the runner gives.
    ln -f -s ~{in_reference_file} reference.fa
    ln -f -s ~{in_reference_index_file} reference.fa.fai
      
    bamleftalign \
        < ~{in_bam_file} \
        > ~{out_prefix}.left_shifted.bam \
        --fasta-reference reference.fa \
        --compressed
    #samtools index -b ~{out_prefix}.left_shifted.bam ~{out_prefix}.left_shifted.bam.bai
  >>>

  output {
    File output_bam_file = "~{out_prefix}.left_shifted.bam"
    #File output_bam_index_file = "~{out_prefix}.left_shifted.bam.bai"
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks : runenv.disks
  }
}
