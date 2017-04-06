#!/bin/bash

# Make some more custom UCSC tracks for Smith, using the new Track Generator script
# https://github.com/stevekm/UCSC-custom-track-generator

# BigWigs
bw_url="https://server.com/ChIP-Seq/alignment/bigwigs/"
bigwig_dir="/data/sequencing/ChIP-Seq/alignment/bigwigs"
D_bw_params="bigwig-green.txt"
R_bw_params="bigwig-purple.txt"
D_out="D_bigwig_tracks.txt"
R_out="R_bigwig_tracks.txt"


(
find "$bigwig_dir" -name '*.bw' -name '*-D-H3K27AC*' | sort | xargs ./make-tracks.py -url "$bw_url" -p "$D_bw_params" -s # -o "$D_out"
find "$bigwig_dir" -name '*.bw' -name '*-R-H3K27AC*' | sort | xargs ./make-tracks.py -url "$bw_url" -p "$R_bw_params" -s # -o "$R_out"
) | sort > Smith_BigWig_tracks.txt

# Super Enhancers
SE_dir="/data/sequencing/ChIP-Seq/superenhancer/ROSE"
SE_url="https://server.com/ChIP-Seq/superenhancer/ROSE/"
SE_params="bigBed-grey.txt"
SE_out="Smith_SE_bigBed_tracks.txt"
find "$SE_dir" -name "*.bigbed" | sort | xargs ./make-tracks.py -url "$SE_url" -p "$SE_params" -o "$SE_out"


# DiffBind Peaks Up Down per Mark
# # all the peaks are the same for every patient so just use these
hg19_chrom_sizes="/ifs/home/kellys04/software/hg19.chrom.sizes"
DiffBed_external_dir="/data/sequencing/ChIP-Seq/diffbind-peaks/per-mark"
DiffBed_external_bed_dir="${DiffBed_external_dir}/bed"
mkdir -p "$DiffBed_external_bed_dir"

DiffBind_Up="/ifs/home/kellys04/projects/ChIpSeq_2016-12-31/project_notes/integrated_analysis/output_per-mark_geneExprTOP-5000_diffRatio-1.5x_methylRatio-1.5x/ABC/H3K27AC_Up.bed"
DiffBind_Up_out_bed="${DiffBed_external_bed_dir}/$(basename "$DiffBind_Up")"
bigBed_name=$(basename "$DiffBind_Up")
bigBed_name="${bigBed_name%%.bed}.bigbed"
DiffBind_Up_out_bigbed="${DiffBed_external_dir}/${bigBed_name}"
/bin/cp -v "$DiffBind_Up" "$DiffBind_Up_out_bed"
bedToBigBed "$DiffBind_Up" "$hg19_chrom_sizes" "${DiffBind_Up_out_bigbed}"

DiffBind_Down="/ifs/home/kellys04/projects/ChIpSeq_2016-12-31/project_notes/integrated_analysis/output_per-mark_geneExprTOP-5000_diffRatio-1.5x_methylRatio-1.5x/ABC/H3K27AC_Down.bed"
DiffBind_Down_out_bed="${DiffBed_external_bed_dir}/$(basename "$DiffBind_Down")"
bigBed_name=$(basename "$DiffBind_Down")
bigBed_name="${bigBed_name%%.bed}.bigbed"
DiffBind_Down_out_bigbed="${DiffBed_external_dir}/${bigBed_name}"
/bin/cp -v "$DiffBind_Down" "$DiffBind_Down_out_bed"
bedToBigBed "$DiffBind_Down" "$hg19_chrom_sizes" "${DiffBind_Down_out_bigbed}"


# DiffBind Peaks per Patient
# need to copy and convert to bigBed format
hg19_chrom_sizes="/ifs/home/kellys04/software/hg19.chrom.sizes"
DiffBed_source_dir="/ifs/home/kellys04/projects/ChIpSeq_2016-12-31/project_notes/integrated_analysis/output_per-patient"
DiffBed_external_dir="/data/sequencing/ChIP-Seq/diffbind-peaks/per-patient"
mkdir -p "$DiffBed_external_dir"


find "${DiffBed_source_dir}" -path "*1.5x_diffRatio-1.5x_methylRatio-1.5x*" -path "*H3K27AC*" -name "Diffbind_Up.bed" | while read item; do
    # /ifs/home/kellys04/projects/ChIpSeq_2016-12-31/project_notes/integrated_analysis/output_per-patient/output_per-patient_geneExprUP-1.5x_diffRatio-1.5x_methylRatio-1.5x/ABC/H3K27AC/Diffbind_Up.bed
    (
    # copy the BED file over to the external dir
    sample_ID="$(basename $(dirname $(dirname "$item") ) )"
    mark_ID="$(basename $(dirname "$item") )"
    sample_output_name="${sample_ID}-${mark_ID}-$(basename "$item")"
    sample_external_bed_dir="${DiffBed_external_dir}/bed"
    mkdir -p "$sample_external_bed_dir"
    output_bed="${sample_external_bed_dir}/${sample_output_name}"
    /bin/cp -v "$item" "$output_bed"

    # convert it to bigBed format
    sample_bigBed="${sample_output_name%%.bed}.bigbed"
    output_bigBed="${DiffBed_external_dir}/${sample_bigBed}"

    bedToBigBed "$item" "$hg19_chrom_sizes" "${output_bigBed}"
    )
done

find "${DiffBed_source_dir}" -path "*1.5x_diffRatio-1.5x_methylRatio-1.5x*" -path "*H3K27AC*" -name "Diffbind_Down.bed" | while read item; do
    # /ifs/home/kellys04/projects/ChIpSeq_2016-12-31/project_notes/integrated_analysis/output_per-patient/output_per-patient_geneExprUP-1.5x_diffRatio-1.5x_methylRatio-1.5x/ABC/H3K27AC/Diffbind_Up.bed
    (
    # copy the BED file over to the external dir
    sample_ID="$(basename $(dirname $(dirname "$item") ) )"
    mark_ID="$(basename $(dirname "$item") )"
    sample_output_name="${sample_ID}-${mark_ID}-$(basename "$item")"
    sample_external_bed_dir="${DiffBed_external_dir}/bed"
    mkdir -p "$sample_external_bed_dir"
    output_bed="${sample_external_bed_dir}/${sample_output_name}"
    /bin/cp -v "$item" "$output_bed"

    # convert it to bigBed format
    sample_bigBed="${sample_output_name%%.bed}.bigbed"
    output_bigBed="${DiffBed_external_dir}/${sample_bigBed}"

    bedToBigBed "$item" "$hg19_chrom_sizes" "${output_bigBed}"
    )
done

# DiffBind bigBed tracks
bigbed_url="https://server.com/ChIP-Seq/diffbind-peaks/per-patient/"
bigbed_dir="/data/sequencing/ChIP-Seq/diffbind-peaks/per-patient"
DiffBind_up_params="bigBed-red.txt"
DiffBind_down_params="bigBed-blue.txt"

(
find "$bigbed_dir" -name "*Diffbind_Down.bigbed" | xargs ./make-tracks.py -url "$bigbed_url" -p "$DiffBind_down_params" -s
find "$bigbed_dir" -name "*Diffbind_Up.bigbed" | xargs ./make-tracks.py -url "$bigbed_url" -p "$DiffBind_up_params" -s
) | sort > Smith_DiffBind_bigBed_tracks.txt

# followed by a lot of manual editting of the track groupings in Excel
