version development

import "../../structs/runenv.wdl"

task run_realigner_target_creator {
  input {
    File in_bam_file
    File in_bam_index_file
    File in_reference_file
    File in_reference_index_file
    File in_reference_dict_file
    Int in_expansion_bases
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

    CONTIG_ID=`head -1 < <(samtools view input_bam_file.bam) | cut -f3`
        
    # Reference and its index must be adjacent and not at arbitrary paths
    # the runner gives.
    ln -f -s "~{in_reference_file}" reference.fa
    ln -f -s "~{in_reference_index_file}" reference.fa.fai
    # And the dict must be adjacent to both
    ln -f -s "~{in_reference_dict_file}" reference.dict

    java -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator \
      --remove_program_records \
      -drf DuplicateRead \
      --disable_bam_indexing \
      -nt ~{runenv.cpu} \
      -R reference.fa \
      -L ${CONTIG_ID} \
      -I input_bam_file.bam \
      --out forIndelRealigner.intervals

      awk -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' forIndelRealigner.intervals > ~{out_prefix}.intervals.bed

      if [ ~{in_expansion_bases} -gt 0 ]; then
        bedtools slop -i ~{out_prefix}.intervals.bed -g "~{in_reference_index_file}" -b "~{in_expansion_bases}" > ~{out_prefix}.intervals.widened.bed
        mv ~{out_prefix}.intervals.widened.bed ~{out_prefix}.intervals.bed
      fi
    >>>

    output {
      File output_target_bed_file = "~{out_prefix}.intervals.bed"
    }

    runtime {
      docker: runenv.docker
      cpu: runenv.cpu
      memory: "~{runenv.memory} GB"
    }
}
