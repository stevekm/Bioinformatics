#!/bin/bash

# submit the 'ccle_make_signatures.R' script to qsub on the HPC to make them all in parrallel since there are a lot of them


# ~~~~~ SETUP ~~~~~ #
ccle_file="ccle2maf_081117.txt"
signatures_script="ccle_make_signatures.R"

qsub_logdir="logs"
mkdir -p "$qsub_logdir"

signatures_dir="signatures"
mkdir -p "$signatures_dir"


sleep_counter=0 # counter to track sleep_count
sleep_limit=200 # the value to sleep on



# ~~~~~ FUNCTIONS ~~~~~ #
run_script () {
    # submit the script to the HPC cluster
    local sampleID="$1"
    local job_name="job_${sampleID}"
    echo "$job_name"
    
    
     qsub -wd $PWD -o :${qsub_logdir}/ -e :${qsub_logdir}/ -j y -N "$job_name" <<E0F
module unload r
# module load r/3.2.3
module load r/3.3.0

set -x
Rscript "${signatures_script}" "${sampleID}"



E0F
    
}


sleep_count () {
    # the count to sleep on
    local sleep_limit="$sleep_limit"

    if [ "$sleep_counter" -eq "$sleep_limit" ] ; then
        printf "\n\nLimit reached. Sleeping....\n\n" 
        sleep 5
        sleep_counter=0
    else
        sleep_counter=$((sleep_counter+1))
    fi
}

# ~~~~~ RUN ~~~~~ #

# iterate over the unique sample IDs in the file
tail -n +2 "$ccle_file" | cut -f16 | sort -u | while read sampleID; do
    # skip empty strings
    if [ ! -z "$sampleID" ] ; then
        run_script "$sampleID"
        sleep_count
    fi
done