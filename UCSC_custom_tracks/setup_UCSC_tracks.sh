#!/bin/bash

# a workflow for setting up some custom UCSC tracks for some ChIP-Seq peaks

external_bed_dir="/ifs/data/sequence/ChIP-Seq/bed"
external_bigbed_dir="/ifs/data/sequence/ChIP-Seq/bigbed"
peaks_dir="/ifs/home/kellys04/projects/ChIpSeq_2016-11-02/pipeline/peaks/results/peaks.by_sample.macs_broad/align.by_sample.bowtie2"
TAD_dir="/ifs/home/kellys04/projects/ChIpSeq_2016-11-02/project_notes/UCSC_browser_sessions/TAD_domains"
hg19_chrom_sizes="/ifs/home/kellys04/software/hg19.chrom.sizes"
bed_to_bigbed="/ifs/home/kellys04/software/UCSC/bedToBigBed"

sort_bed () {
    local bed_file="$1"
    local tmp_file="$(dirname "$bed_file")/tmp"
    sort -k1,1 -k2,2n "$bed_file" > "$tmp_file" && /bin/mv "$tmp_file" "$bed_file"
}

clean_bed () {
    local bed_file="$1"
    local tmp_file="$(dirname "$bed_file")/tmp"
    cat "$bed_file" | cut -f1-3 > "$tmp_file" && /bin/mv "$tmp_file" "$bed_file"
}


find "$peaks_dir" -name "peaks.bed" -path "*CTCF*" | while read item; do
    sampleID="$(basename $(dirname "$item"))"
    output_file="${external_bed_dir}/${sampleID}.bed"
    /bin/cp -v "$item" "$output_file"
    clean_bed "$output_file"
done

find "$TAD_dir" -name "*.bed" ! -path "*__MACOSX*" | while read item; do
    sampleID="$(basename $(dirname "$item"))"
    output_file="${external_bed_dir}/${sampleID}.bed"
    /bin/cp -v "$item" "$output_file"
done

# convert BED to BigBED
find "$external_bed_dir" -name "*.bed" | while read item; do
    item_name="$(basename "$item")"
    output_file="${external_bigbed_dir}/${item_name%%.bed}.bigbed"
    echo "$output_file"
    sort_bed "$item"
    $bed_to_bigbed "$item" "$hg19_chrom_sizes" "$output_file"
done


# get the UCSC track script
# git clone https://github.com/NYU-BFX/UCSC-custom-track-generator.git

external_bigbed_url="https://server.edu/ChIP-Seq/bigbed/"
tracks_script="/ifs/home/kellys04/projects/ChIpSeq_2016-11-02/project_notes/UCSC_browser_sessions/UCSC-custom-track-generator/make-tracks.py"
external_bigbed_dir="/ifs/data/sequence/ChIP-Seq/bigbed"

find "$external_bigbed_dir" -name "*.bigbed" | sort | xargs $tracks_script -url "$external_bigbed_url"

# UCSC Tracks saved to file:
# UCSC_custom_tracks-2017-05-11-17-27-31.txt

# load on https://genome.ucsc.edu/goldenpath/help/customTrack.html
