version development

import "../../structs/runenv.wdl"

task run_realigner_target_creator {
  input {
    File bam
    File bai
    File reference_fasta
    File reference_fai
    File reference_dict
    Int expand_bases = 0
    RunEnv runenv
  }

  String out_prefix = basename(bam, ".bam")
  String intervals = "~{out_prefix}.intervals"
  String targets = "~{out_prefix}.targets.bed"
  String expanded_targets = "~{out_prefix}.targets.expanded.bed"
  command <<<
    # Set the exit code of a pipeline to that of the rightmost command 
    # to exit with a non-zero status, or zero if all commands of the pipeline exit 
    set -euo pipefail
    # echo each line of the script to stdout so we can see what is happening 
    set -o xtrace
    #to turn off echo do 'set +o xtrace' 

    ln -f -s ~{bam} input_bam_file.bam
    ln -f -s ~{bai} input_bam_file.bam.bai

    #CONTIG_ID=`head -1 < <(samtools view input_bam_file.bam) | cut -f3`
        
    # Reference and its index must be adjacent and not at arbitrary paths
    # the runner gives. The dict must be adjacent to both
    ln -f -s "~{reference_fasta}" reference.fa
    ln -f -s "~{reference_fai}" reference.fa.fai
    ln -f -s "~{reference_dict}" reference.dict
    ls -alFh /usr/GenomeAnalysisTK.jar

    # GATK 3.5
    java -Xmx~{runenv.memory - 2}G -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator \
      --remove_program_records \
      -drf DuplicateRead \
      --disable_bam_indexing \
      -nt ~{runenv.cpu} \
      -R reference.fa \
      -I input_bam_file.bam \
      --out ~{intervals}
      #-L ${CONTIG_ID} \

      awk -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' ~{intervals} > ~{targets}
      bedtools slop -i "~{targets}" -g reference.fa.fai -b "~{expand_bases}" > "~{expanded_targets}"
    >>>

  output {
    File intervals = intervals
    File targets = targets
    File expanded_targets = expanded_targets
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "~{runenv.memory} GB"
  }
}
