version development

import "../../../structs/runenv.wdl"

task run_verifybamid {
  input {
    String sample
    File bam
    File bai
    File ref_fasta
    File ref_fai
    File resource
    RunEnv runenv
  }

  String svdprefix = basename(resource, ".tar")
  String output_file = "~{sample}.verifybamid2.txt"
  command <<<
    ln ~{bam} ./
    ln ~{bai} ./
    ln ~{ref_fasta} ./
    ln ~{ref_fai} ./
    mkdir ./resource
    tar xvvf ~{resource} -C resource
    regions_file=$(find resource -name '*.regions')
    # PileUp
    echo "Running mpileup..."
    temp_d=$(mktemp -d)
    trap "rm -rf ${temp_d}" EXIT
    split -l 1000 "${regions_file}" "${temp_d}/regions."
    pileup_file="~{sample}.~{svdprefix}.tsv"
    echo "PileUp file: ${pileup_file}"
    for regions_sub_file in $(find "${temp_d}" -name 'regions.*' ); do
        echo "Regions sub file: ${regions_sub_file}"
        xargs -P ~{runenv.cpu} -I {} bash -c "samtools mpileup -B -a -s -q 1 -r {} -f ~{basename(ref_fasta)} -o ${temp_d}/mpileup_{} ~{basename(bam)} 2>/dev/null" < "${regions_sub_file}" || exit 1
        mpileup_files=(${temp_d}/mpileup_*)
        IFS=$'\n' sorted_files=($(sort -V <<<"${mpileup_files[*]}"))
        unset IFS
        for filename in ${sorted_files[@]}; do
            cat "${filename}" >> "${pileup_file}" || exit 1
            rm -f "${filename}"
        done
    done
    echo "PileUp complete!"
    # VBI
    echo "Running VerifyBamID..."
    VerifyBamID \
      --PileupFile "${pileup_file}" \
      --Reference "~{basename(ref_fasta)}" \
      --SVDPrefix "resource/~{svdprefix}" \
      --NumThread "~{runenv.cpu}" | tee "~{output_file}"
    echo "VerifyBamID output file: ~{output_file}"
  >>>

  output {
    File contamination_report = glob("~{output_file}")[0]
    File ancestry = glob("result.Ancestry")[0]
    File selfsm = glob("result.selfSM")[0]
  }

  runtime {
    docker: runenv.docker
    cpu: runenv.cpu
    memory: "${runenv.memory} GB"
    #disks: runenv.disks
  }
}
