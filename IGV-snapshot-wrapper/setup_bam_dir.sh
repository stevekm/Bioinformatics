#!/bin/bash

proj_dir="/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2017-02-11/project_notes/IGV-snapshots"
align_dir="/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2017-02-11/pipeline/align/results"
bam_outdir="${proj_dir}/bam_links"
mkdir -p "$bam_outdir"

find "$align_dir" -type f -name "*.bam" | while read item; do
    #
    sampleID="$(basename $(dirname $item))"
    sample_group="$(basename $(dirname $(dirname $item)))"
    sample_name="$(echo "$sampleID" | cut -d '-' -f1)"

    printf "%s\t%s\n" "$sampleID" "$sample_group"

    sample_outdir="${bam_outdir}/${sample_group}/${sample_name}"
    mkdir -p "$sample_outdir"
    (
    cd "$sample_outdir"
    ln -fs "$item" "${sampleID}.bam"
    )

done


find "$align_dir" -type f -name "*.bam.bai" | while read item; do
    #
    sampleID="$(basename $(dirname $item))"
    sample_group="$(basename $(dirname $(dirname $item)))"
    sample_name="$(echo "$sampleID" | cut -d '-' -f1)"

    printf "%s\t%s\n" "$sampleID" "$sample_group"

    sample_outdir="${bam_outdir}/${sample_group}/${sample_name}"
    mkdir -p "$sample_outdir"
    (
    cd "$sample_outdir"
    ln -fs "$item" "${sampleID}.bam.bai"
    )

done

