# This is for long reads using the pile up method if we cannot get downsampled bams to work

version development

import "../../../tasks/Utility/Utils.wdl"
import "../../../tasks/Utility/BAMutils.wdl" as BU
import "../../../tasks/QC/Contamination.wdl"

# this is a model that other sub-workflows can potentially follow,
# i.e. define a custom struct so that super workflows can use pre-defined JSON files
struct VBID2_config {
    File genotyping_sites

    Boolean is_hgdp_sites
    Boolean is_100k_sites

    Boolean disable_baq

    Int max_retries

    String? tech
}

workflow LongReadsContaminationEstimation {
    meta {
        desciption:
        "Estimate the cross-individual contamination level of a GRCh38 bam."
    }

    input {
        File bam
        File bai
        String tech
        File ref_map_file

        File gt_sites_bed
        Boolean is_hgdp_sites
        Boolean is_100k_sites

        Boolean disable_baq

        String disk_type

        Int max_retries
    }

    parameter_meta {
        # input:
        gt_sites_bed:     "Bed file holding the genotyping sites."
        is_hgdp_sites:    "Provided BED is HGDP genotyping sites."
        is_100k_sites:    "Provided BED is 100k genotyping sites, not 10k sites."
        disable_baq:      "If turned on, BAQ computation will be disabled (faster operation)."

        tech: "technology used for generating the data; accepted value: [ONT, Sequel, Revio]"

        max_retries: "Because of the strange samtools failures reading from NAS storage, we should make multiple attempts to get away from the trasient errors. If after the max retries, we still get those failures, this task will fail."
    }

    # if the coverage is too low, the tool errors out (and the data won't bring much value anyway)
    # here we guard against it by using bam file size, with a emperically-determined threshold
    Map[String, Int] bam_threshold_per_tech = {'ONT': 450, 'Revio': 150, 'Sequel': 250} # this value is technology dependent
    Int bam_file_threshold = bam_threshold_per_tech[tech]

    if (bam_file_threshold > ceil(size(bam, "MiB"))) {
        Float extreme_low_cov_val = 200  # completely arbitrary
    }

    if (bam_file_threshold <= ceil(size(bam, "MiB"))) {
        # quickly change to pileup
        Map[String, String] ref_map = read_map(ref_map_file)

        Int scaleup_factor = 20
        call BU.BamToRelevantPileup as Pileup {
            input:
                bam = bam,
                bai = bai,
                bed = gt_sites_bed,
                ref_fasta = ref_map['fasta'],
                disable_baq = disable_baq,
                disk_type = disk_type,
                max_retries = max_retries
        }

        call Contamination.VerifyBamID {
            input: pileup = Pileup.pileups, ref_fasta = ref_map['fasta'], is_hgdp_sites = is_hgdp_sites, is_100k_sites = is_100k_sites
        }
    }

    output {
        Float contamination_est = select_first([VerifyBamID.contamination_est, extreme_low_cov_val])
    }
}

task BamToRelevantPileup {
    meta {
        desciption:
        "Chop up a GRCh38 BAM by chromosome and further subset into requested genotyping sites; then convert to pileup format. See also task GetPileup."
    }

    parameter_meta {
        bam: {localization_optional: true}
        bed: "sites where pileup is needed"
    }
    input {
        File bam
        File bai
        File bed
        File ref_fasta
        Boolean disable_baq

        String disk_type = "SSD"
    }
    output {
        File pileups = "pileup.mpileup"
    }

    String baq_option = if disable_baq then '-B' else '-E'

    String base = basename(bam)
    String local_bam = "/cromwell_root/~{base}"

    command <<<
        set -euxo pipefail

        time \
        gcloud storage cp ~{bam} ~{local_bam}
        mv ~{bai} "~{local_bam}.bai"

        # generate bed for parallel conversion
        set +e
        for i in `seq 1 22`;
        do
            grep -w "chr${i}" ~{bed} > "chr${i}.bed";
        done
        grep -w "chrX" ~{bed} > "chrX.bed"
        grep -w "chrY" ~{bed} > "chrY.bed"
        set -e
        rm ~{bed}

        # parallel conversion
        cnt=0
        for bed in $(ls chr*.bed | sort -V); do

            if [[ ! -s ${bed} ]] ; then rm "${bed}" && continue; fi

            prefix=$(echo "${bed}" | awk -F '.' '{print $1}')
            samtools view -h \
                --region-file "${bed}" \
                --write-index \
                -o "${prefix}.bam##idx##${prefix}.bam.bai" \
                ~{local_bam} && \
            samtools mpileup \
                ~{baq_option} \
                -s \
                -q 1 \
                -f ~{ref_fasta} \
                -o "${prefix}.mpileup" \
                "${prefix}.bam" &
            cnt=$((cnt + 1))
            if [[ $cnt -eq ~{cores} ]]; then wait; cnt=0; rm -f chr*bam chr*bai; fi
        done
        wait

        rm -f chr*bam chr*bai
        cat *.mpileup > pileup.mpileup
    >>>

    Int cores = 12
    Int memory = 4 + cores
    Int local_ssd_sz = if size(bam, "GiB") > 150 then 750 else 375
    Int pd_sz = 10 + 2 * ceil(size(bam, "GiB"))
    Int disk_size = if "LOCAL" == disk_type then local_ssd_sz else pd_sz

    runtime {
        cpu:            "~{cores}"
        memory:         "~{memory} GiB"
        disks:          "local-disk ~{disk_size} ~{disk_type}"
        preemptible:    1
        maxRetries:     1
        docker:         "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
}

task VerifyBamID {
    meta {
        desciption: "Uses VerifyBamID2 for human cross-individual contamination estimation. Assumes GRCh38."
    }

    input {
        File pileup
        File ref_fasta
        Boolean is_hgdp_sites
        Boolean is_100k_sites
    }

    String a = if is_hgdp_sites then 'hgdp' else '1000g.phase3'
    String b = if is_100k_sites then '100k' else  '10k'
    String resource_prefix = '~{a}.~{b}.b38.vcf.gz.dat'

    command <<<
        set -eux

        export VERIFY_BAM_ID_HOME='/VerifyBamID'

        time \
        ${VERIFY_BAM_ID_HOME}/bin/VerifyBamID \
            --SVDPrefix ${VERIFY_BAM_ID_HOME}/resource/~{resource_prefix} \
            --Reference ~{ref_fasta} \
            --PileupFile ~{pileup} \
            --NumThread 4 \
        > vbid2.out \
        2> vbid2.log

        cat vbid2.out
        tail -1 vbid2.out | awk -F ':' '{print $2}' | awk '{$1=$1};1' > "est_contam.txt"
    >>>

    output {
        File vbid2_log = "vbid2.log"
        File vbid2_out = "vbid2.out"
        Float contamination_est = read_float("est_contam.txt")
    }

    Int disk_size = 10 + ceil(size(pileup, "GiB"))
    runtime {
        cpu: 4
        memory: "8 GiB"
        disks: "local-disk ~{disk_size} SSD"
        docker: "us.gcr.io/broad-dsp-lrma/verifybamid2:v2.0.1"
    }
}
