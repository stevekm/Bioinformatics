#!/bin/bash

# run the script for making ChIP-Seq profile plots
# submit the R script to the cluster for every sample cause its real slow


# ~~~~~ LOCATIONS ~~~~~ #
output_dir='/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2017-12-31/project_notes/profile_plots/output'
peaks_dir='/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2017-12-31/pipeline/peaks/results/peaks.by_group.macs_broad_nolambda/align.by_sample.bowtie2'
R_script="/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2017-12-31/project_notes/profile_plots/ChIPseeker_profile_plots.R"

# ~~~~~ FUNCTIONS ~~~~~ #
submit_job () {
    local filepath="$1"
    local sampleID="$(basename "$(dirname "$filepath")" )"
    local output_path="${output_dir}/${sampleID}_profile.pdf"
    local qsub_logdir="${output_dir}/logs"
    local job_name="$sampleID"

    mkdir -p "$qsub_logdir"

    # printf '%s %s %s\n' "$sampleID" "$filepath" "$output_path"

    script_command="$R_script $sampleID $filepath $output_path"
    # /ifs/home/kellys04/projects/SmithLab_ChIpSeq_2017-12-31/project_notes/profile_plots/ChIPseeker_profile_plots.R Sample1-R-H3K27AC /ifs/home/kellys04/projects/SmithLab_ChIpSeq_2017-12-31/pipeline/peaks/results/peaks.by_group.macs_broad_nolambda/align.by_sample.bowtie2/Sample1-R-H3K27AC/peaks.bed /ifs/home/kellys04/projects/SmithLab_ChIpSeq_2017-12-31/project_notes/profile_plots/output/Sample1-R-H3K27AC_profile.pdf

    qsub -wd $PWD -o :${qsub_logdir}/ -e :${qsub_logdir}/ -j y -N "$job_name" <<E0F
    set -x
    date +"%Y-%m-%d-%H-%M-%S"
    start=\$(date +%s)

    $script_command

    date +"%Y-%m-%d-%H-%M-%S"
    echo "Duration: \$(((\$(date +%s)-\$start)/60)) minutes"
    set +x
E0F

}


# ~~~~~ RUN ~~~~~ #
find "$peaks_dir" -type f -name "peaks.bed" -print0 | while read -d $'\0' item; do
    submit_job "$item"
done
