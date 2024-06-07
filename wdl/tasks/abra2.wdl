version development

import "../structs/runenv.wdl"

task run_realigner {
  input {
    File in_bam_file
    File in_bam_index_file
    File in_target_bed_file
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

    ln -f -s ~{in_bam_file} input_bam_file.bam
    ln -f -s ~{in_bam_index_file} input_bam_file.bam.bai

    # Reference and its index must be adjacent and not at arbitrary paths
    # the runner gives.
    ln -f -s ~{in_reference_file} reference.fa
    ln -f -s ~{in_reference_index_file} reference.fa.fai

    java -Xmx~{runenv.memory - 2}G -jar /apps/abra2/abra2.jar \
      --targets ~{in_target_bed_file} \
      --in input_bam_file.bam \
      --out ~{out_prefix}.indel_realigned.bam \
      --ref reference.fa \
      --index \
      --threads ~{runenv.cpu}
  >>>

  output {
    File indel_realigned_bam = "~{out_prefix}.indel_realigned.bam"
    File indel_realigned_bam_index = "~{out_prefix}.indel_realigned.bai"
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: runenv.memory + " GB"
    disks : runenv.disks
  }
}
