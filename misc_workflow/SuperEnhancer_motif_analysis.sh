#!/bin/bash

# Motif analysis workflow for Super Enhancer BED file inputs 
# using HOMER and MEME-CHIP, also requires bedtools

project_dir="/ifs/home/kellys04/ChIPSeq_2016-12-31/motifs-SuperEnhancers"
SE_dir="/ifs/data/ChIP-Seq/superenhancer/ROSE/bed/"

output_HOMER="${project_dir}/HOMER"
mkdir -p "$output_HOMER"

output_MEME="${project_dir}/MEME-CHIP"
mkdir -p "$output_MEME"

run_HOMER_motif () {
    # http://homer.salk.edu/homer/ngs/peakMotifs.html
    local input_bed="$1" # ABC-D-H3K27AC_ROSE_superenhancer.bed
    local genome="hg19"
    local motif_options="-size given"
    local sample_ID="$(basename "$input_bed" | cut -d '-' -f1)"
    local sample_status="$(basename "$input_bed" | cut -d '-' -f2)"
    local outdir="${output_HOMER}/${sample_ID}/${sample_status}"
    mkdir -p "$outdir"
    local preparsed_dir="${outdir}/preparsed"
    mkdir -p "$preparsed_dir"
    local logdir="${outdir}/logs"
    mkdir -p "$logdir"
    local job_name="homer-${sample_ID}"
    echo "$sample_ID $sample_status"
    echo "$outdir"
    qsub -wd $outdir -o :${logdir}/ -e :${logdir}/ -j y -pe threaded 4 -N "$job_name" <<E0F
module load homer/v4.6 # set environment
set -x
threads=\$NSLOTS # THREADS=${NSLOTS:=4}
findMotifsGenome.pl $input_bed $genome $outdir/ -p $threads -preparse -preparsedDir ${preparsed_dir}/ -dumpFasta $motif_options # && rm -rf $preparsed_dir
set +x
E0F
}

bed2fasta () {
    # make a fasta from a BED for motif analysis
    local input_bed="$1"
    local sample_ID="$(basename "$input_bed" | cut -d '-' -f1)"
    local sample_status="$(basename "$input_bed" | cut -d '-' -f2)"
    local outdir="${output_MEME}/${sample_ID}/${sample_status}"
    mkdir -p "$outdir"
    local input_bed_basename="$(basename "$input_bed")"
    local output_fasta="${outdir}/${input_bed_basename%%.bed}.fasta"
    local genome_fasta="/local/data/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
    echo "$sample_ID $sample_status $output_fasta"
    (
    module unload bedtools
    module unload gcc
    module load bedtools/2.22.0
    bedtools getfasta -s -fi "$genome_fasta" -bed "$input_bed" -fo "$output_fasta"
    )
}

run_MEME_CHIP () {
    local input_fasta="$1" # ....motifs-SuperEnhancers/MEME-CHIP/ABC/R/ABC-R-H3K27AC_ROSE_superenhancer.fasta
    local sample_ID="$(basename "$input_fasta" | cut -d '-' -f1)"
    local sample_status="$(basename "$input_fasta" | cut -d '-' -f2)"
    local outdir="${output_MEME}/${sample_ID}/${sample_status}"
    mkdir -p "$outdir"
    local input_bed_basename="$(basename "$input_fasta")"
    local JasparDB="/ifs/home/kellys04/data/motifs/motif_databases/JASPAR/JASPAR_CORE_2016_vertebrates.meme"
    local JOLMAdb="/ifs/home/kellys04/data/motifs/motif_databases/EUKARYOTE/jolma2013.meme"
    local MIRBASE_human_db="/ifs/home/kellys04/data/motifs/motif_databases/MIRBASE/Homo_sapiens_hsa.dna_encoded.meme"
    local HOCOMOCO_db="/ifs/home/kellys04/data/motifs/motif_databases/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme"
    local logdir="${outdir}/logs"
    mkdir -p "$logdir"
    local job_name="MEME-CHIP-${sample_ID}"
    echo "$sample_ID $sample_status"
    echo "$outdir"
    qsub -wd $outdir -o :${logdir}/ -e :${logdir}/ -j y -pe threaded 4 -N "$job_name" <<E0F


# # module for this version are currently broken but this is the code I used previously
# module unload gcc
# module unload perl
# module unload python
# module load meme/4.11.1
# meme-chip -meme-p 4 -oc "$outdir" -index-name meme-chip.html -time 300 -order 1 -norand -db "$JOLMAdb" -db "$JasparDB" -db "$MIRBASE_human_db" -db "$HOCOMOCO_db" -meme-mod zoops -meme-minw 6 -meme-maxw 30 -meme-nmotifs 5 -dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 "$input_fasta"

module load meme/4.9.1
set -x
meme-chip --version
meme-chip -meme-p 4 -oc "$outdir" -index-name meme-chip.html -time 300 -db "$JOLMAdb" -db "$JasparDB" -db "$MIRBASE_human_db" -db "$HOCOMOCO_db" -meme-mod zoops -meme-minw 6 -meme-maxw 30 -meme-nmotifs 5 -dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 "$input_fasta"

E0F
}

# HOMER
find "$SE_dir" -name "*.bed" -exec readlink -f {} \; | while read item; do
    run_HOMER_motif "$item"
done

# .bed -> .fasta
find "$SE_dir" -name "*.bed" -exec readlink -f {} \; | while read item; do
    bed2fasta "$item"
done

# MEME-CHiP
find "$output_MEME" -name "*.fasta" -exec readlink -f {} \; | while read item; do
    run_MEME_CHIP "$item"
done
