#!/bin/bash

# find bam files in the dir and convert them to bigwig

conv_bam_to_bigwig () {
    local bamfile="$1"
    local bigwigfile="${bamfile%%.bam}.bw"
    local genome="hg19"
    (
    module load samtools/0.1.19
    module load bedtools/2.22.0
    echo "$bamfile $bigwigfile"
    python -c 'from pybedtools.contrib.bigwig import bam_to_bigwig; import sys; bam_to_bigwig(genome=str(sys.argv[1]), bam=str(sys.argv[2]), output=str(sys.argv[3]))' "$genome" "$bamfile" "$bigwigfile"
    )
}

cd /ifs/home/kellys04/projects/ChIpSeq_2017-07-19/project_notes/the_search_for_peaks/data_links
find . -maxdepth 1 -name "*.bam" | while read item; do
    conv_bam_to_bigwig "$item"
done
