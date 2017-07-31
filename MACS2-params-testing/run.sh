#!/bin/bash

# run the MACS2 settings with the given params settings


# ~~~~~ SETUP ~~~~~ #
proj_dir="/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2017-07-19/project_notes/test-macs2-params"
cd "$proj_dir"

samplesheet="selected-samples.tsv"
params_file="macs2_params.tsv"

output_dir="${proj_dir}/output"
mkdir -p "$output_dir"

log_dir="${proj_dir}/logs"
mkdir -p "$log_dir"

# ~~~~~ FUNCTIONS ~~~~~ #
run_samplesheet () {
    # run entries from each line of the samplesheet
    local samplesheet="$samplesheet"
    cat "$samplesheet" | while read line; do
        if [ ! -z "$line" ]; then
            sampleID="$(echo "$line" | cut -f1)"
            sample_bam="$(echo "$line" | cut -f2)"
            control_bam="$(echo "$line" | cut -f3)"
            if [ ! -z "$sampleID" ] && [ ! -z "$sample_bam" ] && [ ! -z "$control_bam" ]; then
                sample_outdir="${output_dir}/${sampleID}"
                mkdir -p "$sample_outdir"
                run_macs_params "$sampleID" "$sample_bam" "$control_bam" "$sample_outdir"
            fi
        fi
    done
}

run_macs_params () {
    # run MACS2 based on the params on each line of the params sheet
    local sampleID="$1"
    local sample_bam="$2"
    local control_bam="$3"
    local sample_outdir="$4"
    local params_file="$params_file"

    cat "$params_file" | while read line; do
        if [ ! -z "$line" ]; then
            params_name="$(echo "$line" | cut -f1)"
            params="$(echo "$line" | cut -f2)"
            run_macs "$sampleID" "$sample_bam" "$control_bam" "$params_name" "$params" "$sample_outdir"
        fi
    done
}

run_macs () {
    # run MACS2
    local sampleID="$1"
    local sample_bam="$2"
    local control_bam="$3"
    local name="$4"
    local params="$5"
    local sample_outdir="$6"
    local sample_params_name="${sampleID}_${name}"
    # printf "%s %s %s %s %s %s\n" "$sampleID" "$sample_bam" "$control_bam" "$params" "$name" "$params"
    # macs2 callpeak -t inpdirs/align/results/align.by_sample.bowtie2/BVI-D-H3K27AC/alignments.bam -c inpdirs/align/results/align.by_sample.bowtie2/BVI-R-INPUT/alignments.bam --outdir=results/peaks.by_sample.macs_broad/align.by_sample.bowtie2/BVI-D-H3K27AC --name=macs --broad --nomodel --extsize=200 -g hs

    # check if there are already output files
    output_files="$(find "$sample_outdir" -name "*$sample_params_name*" )"

    # Only run if none were found
    if [ -z "${output_files}" ]; then
        (
        qsub -wd $PWD -o :${log_dir}/ -e :${log_dir}/ -j y -N "$sample_params_name" <<E0F
        source venv/bin/activate
        set -x
        macs2 --version
        macs2 callpeak -t "$sample_bam" -c "$control_bam" --outdir="$sample_outdir" --name="$sample_params_name" $params
E0F
    sleep 2 # to avoid file locks, etc.
    )
    fi
}

# ~~~~~ RUN ~~~~~ #
run_samplesheet "$samplesheet"


exit
