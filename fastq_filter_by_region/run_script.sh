#!/bin/bash

# run the scrip on the HPC cluster


# ~~~~~ SETUP ~~~~~ #
# script to run
script="run_script.py"

# dir for qsub logs
qsub_logdir="logs"
mkdir -p "$qsub_logdir"

# file with sample IDs, paths to .fastq and .bam files
samples_file="sample_files.tsv"



# ~~~~~ FUNCTIONS ~~~~~ #

run_script () {
    local sampleID="$1"
    local fastqs="$2"
    local bams="$3"
    local job_name="job_$sampleID"

    cmd="./${script} ${fastqs} ${bams}"

    echo "$cmd"

    set -x

    qsub -wd $PWD -o :${qsub_logdir}/ -e :${qsub_logdir}/ -j y -N "$job_name" <<E0F

module unload python
module load python/2.7

set -x

$cmd

E0F
set +x
echo ""
}



# ~~~~~ RUN ~~~~~ #
tail -n+2 "$samples_file" | while read line; do
    if [ ! -z "$line" ]; then
        sampleID="$(echo "$line" | cut -f1)"
        bams="$(echo "$line" | cut -f2)"
        fastqs="$(echo "$line" | cut -f3)"

        run_script "$sampleID" "$fastqs" "$bams"
    fi
done
