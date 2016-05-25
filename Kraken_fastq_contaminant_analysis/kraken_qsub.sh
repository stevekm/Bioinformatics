#!/bin/bash
set -x
##
## USAGE: kraken_qsub.sh /path/to/fastqc_outdir /path/to/input_file <sampleID>
## This script will run the Kraken fastq contaminant analysis to see what the reads might have come from
## 

# ~~~~~~ permissions ~~~~~~ #
# To make stuff group-writeable (this is what I add to all my shared scripts):
umask 007


# ~~~~~~ script args processing ~~~~~~ #
# # if not enough args, output USAGE info lines (which start with '##') and exit
if (($# != 3)); then
  grep '^##' $0
  exit
fi

# get the arguments
OUTDIR="$1" # outdir
INPUTFILE="$2" # input file
SAMPLEID="$3"
THREADS=${NSLOTS:=6}
echo "Working dir is $PWD"
echo "OUTDIR is $OUTDIR"
echo "INPUTFILE is $INPUTFILE"
echo "SAMPLEID is $SAMPLEID"
echo "THREADS is $THREADS"

# kraken settings
export KRAKEN_DEFAULT_DB="$HOME/ref/Kraken/nt-48G"
export KRAKEN_NUM_THREADS="$THREADS"


# ~~~~~~ Script Logging ~~~~~~~ #
# because you can never have too many logs; make a hard copy of the current script
# get the script file 
zz=$(basename $0)
# get the script dir
za=$(dirname $0)
# for regular script running (no qsub)
# LOG_FILE=log.$(basename $0).$(date -u +%Y%m%dt%H%M%S).${HOSTNAME:-$(hostname)}.$$.${RANDOM}
# for qsub script submission
LOG_FILE=scriptlog.$JOB_NAME.$(date -u +%Y%m%dt%H%M%S).${HOSTNAME:-$(hostname)}.$$.${RANDOM}
echo "This is the log file for the script." > $LOG_FILE
echo -e "\n pwd is:\n$PWD\n" >> $LOG_FILE
# echo -e "\nScript file is:\n$0\n" >> $LOG_FILE # for regular script usage (no qsub)
echo -e "\nScript file is:\n$JOB_NAME\n" >> $LOG_FILE # for qsub
echo -e "\nScript file contents:\n" >> $LOG_FILE
cat "$0" >> $LOG_FILE
# ~~~~~~~~~~~~ #

# ~~~~~~ run command ~~~~~~ #
# Thanks to igordot for this; https://github.com/igordot
zcat $INPUTFILE | head -4000000 | \
$HOME/software/bin/kraken --fastq-input /dev/fd/0 | \
$HOME/software/bin/kraken-report --show-zeros | awk -F $'\t' '$1>0.1' > kraken_contaminant_analysis.${SAMPLEID}.txt

# filter for top hits
cat kraken_contaminant_analysis.${SAMPLEID}.txt | awk -F $'\t' '$1>1.0 && $4!="-"' | cut -f 1,2,4,5,6 > tophit.txt

# version
$HOME/software/bin/kraken --version
