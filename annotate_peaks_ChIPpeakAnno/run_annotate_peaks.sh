#!/bin/bash

# this script will find the peaks files in the input directory
# and run the annotation script on them

project_dir="/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2017-02-11_2/project_notes/peak_annotation_stats"
input_dir="/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2017-02-11_2/pipeline/peaks/results"
peaks_file_basename="peaks.bed"

output_dir="${project_dir}/annotated_peaks"
mkdir -p "$output_dir"

annotation_script="${project_dir}/annotate_peaks.R"

# find all the peaks files to be annotated
find "$input_dir" -name "$peaks_file_basename" | while read item; do
set -x
echo "$item"

input_name="$(echo "$item" | sed -e "s|${input_dir}/||g")"
sample_ID="$(basename $(dirname "$input_name"))"
sample_path="$(dirname "$input_name")"

output_path="${output_dir}/${sample_path}"
log_dir="${output_path}/logs"
mkdir -p "$log_dir"

output_filepath="${output_path}/${sample_ID}_annotated_peaks.tsv"
echo "$output_filepath"

echo ""
qsub -b y -wd "$project_dir" -o :${log_dir}/ -e :${log_dir}/ -N "$sample_ID" Rscript "$annotation_script" "$item" "$output_filepath" # -pe threaded "$job_threads" -l mem_free="$job_mem" -l h_vmem="$job_mem" -l mem_token="$job_mem"

# Rscript "$annotation_script" "$item" "$output_filepath"

set +x
done
