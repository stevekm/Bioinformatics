#!/bin/bash
# print the script commands to stderr on execution; stderr gets saved in the qsub log files
set -x
##
## USAGE: FastQC_Kraken_Krona.sh /path/to/outdir /path/to/input_file.fastq.gz 
## 

# ~~~~~~ permissions ~~~~~~ #
# To make stuff group-writeable
umask 007

# ~~~~~~ script args processing ~~~~~~ #
# # if not enough args, output USAGE info lines (which start with '##') and exit
if (($# != 2)); then
  grep '^##' $0
  exit
fi

# ~~~~~~ params ~~~~~~ #
# these need to be installed on the cluster
export KRAKEN_DEFAULT_DB="$HOME/ref/Kraken/nt-48G"
export KRAKEN_NUM_THREADS="$THREADS"
module load fastqc/0.11.4 # use this if FastQC is not in your $PATH

# get the script arguments
OUTDIR="$1" # outdir
mkdir -p "$OUTDIR"
FASTQ_R1="$2" # Read 1 input file


# THREADS is defined by qsub; if not, set with 8 threads
THREADS=${NSLOTS:=2} 
echo "OUTDIR is $OUTDIR"
echo "FASTQ_R1 is $FASTQ_R1"
echo "THREADS is $THREADS"


# ~~~~~~ # ~~~~~~  QUALITY CONTROL # ~~~~~~ # ~~~~~~ # ~~~~~~ 
# FastQC; check raw data quality
fastqc --threads "$THREADS" --nogroup --outdir "$OUTDIR" "$FASTQ_R1"


# Kraken; check for potential contaminants
output_file="$(basename $FASTQ_R1)"
# # only use the first 1,000,000 reads !! # 4,000,000
zcat $FASTQ_R1 | head -400 | \
$HOME/software/bin/kraken --fastq-input /dev/fd/0 | \
tee >(cut -f2,3 > "${OUTDIR}/${output_file}_krona.input") | \
$HOME/software/bin/kraken-report --show-zeros | awk -F $'\t' '$1>0.1' > "${OUTDIR}/${output_file}_kraken_contaminant_analysis.txt"

# filter for top hits
cat "${OUTDIR}/${output_file}_kraken_contaminant_analysis.txt" | awk -F $'\t' '$1>1.0 && $4!="-"' | cut -f 1,2,4,5,6 > "${OUTDIR}/${output_file}_kraken_top_hits.txt"

# get the software version
# $HOME/software/bin/kraken --version
# Kraken version 0.10.6-unreleased
# Copyright 2013-2015, Derrick Wood (dwood@cs.jhu.edu)

# create Krona plot
~/software/bin/Krona/bin/ktImportTaxonomy "${OUTDIR}/${output_file}_krona.input" -o "${OUTDIR}/${output_file}_krona_plot.html"

# use the editted raster script with PhantomJS 
~/software/bin/phantomjs ~/software/bin/rasterize_1000.js "${OUTDIR}/${output_file}_krona_plot.html" "${OUTDIR}/${output_file}_krona_snap.png"
# https://github.com/ariya/phantomjs
# http://phantomjs.org/download.html 