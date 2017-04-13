#!/bin/bash

# Motif analysis workflow for peaks that overlap the Super Enhancer regions

project_dir="/ifs/home/kellys04/projects/ChIpSeq_2016-12-31/project_notes/motifs-SuperEnhancers"
SE_dir="/ifs/data/sequence//ChIP-Seq/superenhancer/ROSE/bed/"
peak_dir="/ifs/home/kellys04/projects/ChIpSeq_2016-12-31/pipeline/peaks/results/peaks.by_sample.macs_broad/align.by_sample.bowtie2"
output_dir="${project_dir}/Peak_SE_overlap_motifs"
mkdir -p "$output_dir"


find_peak_overlaps () {
    local peaks_bed="$1"
    local SE_bed="$2"
    local output_bed="$3"
    (
    module load bedtools/2.26.0
    bedtools intersect -a "$peaks_bed" -b "$SE_bed" -wa | sort -u > "$output_bed"
    )
}


bed2fasta () {
    # make a fasta from a BED for motif analysis
    local input_bed="$1"
    local output_fasta="$2"
    local genome_fasta="/local/data/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
    (
    module unload bedtools
    module unload gcc
    module load bedtools/2.22.0
    bedtools getfasta -s -fi "$genome_fasta" -bed "$input_bed" -fo "$output_fasta"
    )
}


run_MEME_CHIP () {
    local input_fasta="$1" 
    local outdir="$2"
    local logdir="${outdir}/logs"
    mkdir -p "$logdir"
    local sample_ID="$3"
    local job_name="MEME-${sample_ID}"
    local JasparDB="/ifs/home/kellys04/data/motifs/motif_databases/JASPAR/JASPAR_CORE_2016_vertebrates.meme"
    local JOLMAdb="/ifs/home/kellys04/data/motifs/motif_databases/EUKARYOTE/jolma2013.meme"
    local MIRBASE_human_db="/ifs/home/kellys04/data/motifs/motif_databases/MIRBASE/Homo_sapiens_hsa.dna_encoded.meme"
    local HOCOMOCO_db="/ifs/home/kellys04/data/motifs/motif_databases/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme"
    qsub -wd $outdir -o :${logdir}/ -e :${logdir}/ -j y -pe threaded 4-16 -N "$job_name" <<E0F
module unload python
module load meme/4.11.2

set -x
threads=\$NSLOTS # THREADS=${NSLOTS:=8}
meme-chip --version
meme-chip -meme-p \$threads -oc "$outdir" -index-name meme-chip.html -time 300 -order 1 -norand -db "$JOLMAdb" -db "$JasparDB" -db "$MIRBASE_human_db" -db "$HOCOMOCO_db" -meme-mod zoops -meme-minw 6 -meme-maxw 30 -meme-nmotifs 5 -dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 "$input_fasta"

set +x
E0F
}

# for every Super Enhancer BED file, find the corresponding peaks file, convert to fasta, run MEME-CHIP
find "${SE_dir}" -name "*.bed" -exec readlink -f {} \; | while read sample_SE_bed; do
    sample_ID="$(basename "$sample_SE_bed" | cut -d '_' -f1)"
    echo "$sample_ID"
    sample_outdir="${output_dir}/${sample_ID}"
    mkdir -p "$sample_outdir"
    sample_peaks_bed="$(find "$peak_dir" -path "*${sample_ID}*" -name "peaks.bed" | head -1)"
    echo "$sample_SE_bed"
    sample_peak_SE_overlaps_bed="${sample_outdir}/${sample_ID}_peak_overlaps.bed"
    echo "$sample_peak_SE_overlaps_bed"
    sample_peak_SE_overlaps_fasta="${sample_peak_SE_overlaps_bed%%.bed}.fasta"
    echo "$sample_peak_SE_overlaps_fasta"
    sample_MEME_outdir="${sample_outdir}/MEME"

    find_peak_overlaps "$sample_peaks_bed" "$sample_SE_bed" "$sample_peak_SE_overlaps_bed"
    bed2fasta "$sample_peak_SE_overlaps_bed" "$sample_peak_SE_overlaps_fasta"
    run_MEME_CHIP "$sample_peak_SE_overlaps_fasta" "$sample_MEME_outdir" "$sample_ID"
done
