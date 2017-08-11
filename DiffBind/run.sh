#!/bin/bash

# run the custom diffbind

# ~~~~~ SETUP ~~~~~ #
proj_dir="/ifs/home/kellys04/projects/SmithLab_ChIPSeq_2017-12-31/project_notes/diffbind-custom"
output_parent_dir="${proj_dir}/output"
mkdir -p "$output_parent_dir"
cd "$proj_dir"

diffbind_script="${proj_dir}/pipeline-chipseq-diffbind.R"
# ./pipeline-chipseq-diffbind.R output/ diffbind-sample-sheet.csv hg19 foo

diffbind_samplesheet_base="${proj_dir}/diffbind-sample-sheet-"
diffbind_samplesheet_ext=".csv"
marks="H3K9ME3 H3K9AC H3K4ME3 H3K27ME3 H3K27AC CTCF"



# ~~~~~ FUNCTIONS ~~~~~ #
timestamp () {
    printf "%s" "\$(date +"%Y-%m-%d-%H-%M-%S")"
}

run_diffbind () {
    local mark="$1"
    # /ifs/home/kellys04/projects/SmithLab_ChIPSeq_2017-12-31/project_notes/diffbind-custom/diffbind-sample-sheet-H3K27ME3.csv
    local samplesheet="${diffbind_samplesheet_base}${mark}${diffbind_samplesheet_ext}"

    echo "$samplesheet"

    if [ -f "$samplesheet" ]; then
        # setup output location
        local output_dir="${output_parent_dir}/${mark}"
        local samplesheet_basename="$(basename "$samplesheet")"
        local samplesheet_copy="${output_dir}/${samplesheet_basename}"
        mkdir -p "$output_dir"
        /bin/cp -v "$samplesheet" "$samplesheet_copy"

        # submit qsub job
        qsub -wd $PWD -j y -N "${mark}-diffbind" -o :${output_dir}/ -e :${output_dir}/ <<E0F
set -x
printf "qsub job started at: %s" "$(timestamp)"

# ./pipeline-chipseq-diffbind.R output/ diffbind-sample-sheet.csv hg19 foo
${diffbind_script} "${output_dir}" "$samplesheet_copy" hg19
E0F

sleep 1
    fi
}


# ~~~~~ RUN ~~~~~ #
for mark in $marks; do
    run_diffbind "$mark"
done
