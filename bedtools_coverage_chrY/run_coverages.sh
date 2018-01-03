#!/bin/bash

# create coverages of all chrY bases in the alignments
# http://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html

# ~~~~~ SETUP ~~~~~ #
hg19_chrY="/local/data/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/chrY.fa"
qsub_logdir="${PWD}/logs"
mkdir -p "$qsub_logdir"
output_dir="${PWD}/output"
mkdir -p "$output_dir"


# ~~~~~ FUNCTIONS ~~~~~ #
run_cov () {
    # run bedtools genomecov on a sample
    local sampleID="$1"
    local job_name="job_${sampleID}"
    local bamfile="$2"
    
    # all covereage values for every base
    local output_cov="${output_dir}/${sampleID}.chrY.cov.txt"
    # coverage for all contiguously covered regions
    local output_bedgraph="${output_dir}/${sampleID}.chrY.all.bedgraph"
    # the same as above but without 0 regions
    local output_good_bedgraph="${output_dir}/${sampleID}.chrY.bedgraph"
    
    qsub -wd $PWD -o :${qsub_logdir}/ -e :${qsub_logdir}/ -j y -N "$job_name" <<E0F

module load bedtools/2.26.0
set -x

bedtools genomecov -ibam "$bamfile" -g "$hg19_chrY" -bg | grep -w 'chrY' > "$output_good_bedgraph"
bedtools genomecov -ibam "$bamfile" -g "$hg19_chrY" -bga | grep -w 'chrY' > "$output_bedgraph"
bedtools genomecov -ibam "$bamfile" -g "$hg19_chrY" -d | grep -w 'chrY' > "$output_cov"

E0F

}


# ~~~~~ RUN ~~~~~ #
# find all the alignment files in the output directory
find "$(readlink -f BAM-BWA)" -name "*.bam" | while read file; do
    sampleID="$(basename "$file")"
    sampleID="${sampleID%%.bam}"
    
    printf "%s %s\n\n" "$sampleID" "$file"
    run_cov "$sampleID" "$file"
done
