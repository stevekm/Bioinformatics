#!/bin/bash

##
## USAGE: fastqc_qsub.sh /path/to/fastqc_outdir /path/to/input_file
## This script will run FastQC; this script needs to be submitted with qsub
## 


# ~~~~~~ script args ~~~~~~ #
# # if not enough args, output USAGE info lines (which start with '##') and exit
if (($# != 2)); then
  grep '^##' $0
  exit
fi


# make sure that the correct number of script arguments were used
# if [ $# != 2 ] # if not enough args provided
# then
#   grep '^##' $0 # print out lines from the script that start with '##' and exit
#   exit
# fi


# get the arguments
FASTQC_DIR="$1" # outdir
FASTQ1="$2" # input file
THREADS=$NSLOTS #  from qsub

# set number of threads; check if NSLOTS from qsub is set, otherwise use 4
THREADS=${NSLOTS:=4}

# ~~~~~~ FASTQC Params ~~~~~~ #
module load fastqc/0.11.4
# ~~~~~~ # ~~~~~~ # ~~~~~~ #


# ~~~~~~ Script Logging ~~~~~~~ #
# get the script file 
zz=$(basename $0)
# get the script dir
za=$(dirname $0)
# for regular script running (no qsub)
# LOG_FILE=log.$(basename $0).$(date -u +%Y%m%dt%H%M%S).${HOSTNAME:-$(hostname)}.$$.${RANDOM}
# for qsub script submission
LOG_FILE=log.$JOB_NAME.$(date -u +%Y%m%dt%H%M%S).${HOSTNAME:-$(hostname)}.$$.${RANDOM}
echo "This is the log file for the script." > $LOG_FILE
echo -e "\n pwd is:\n$PWD\n" >> $LOG_FILE
# echo -e "\nScript file is:\n$0\n" >> $LOG_FILE # for regular script usage (no qsub)
echo -e "\nScript file is:\n$JOB_NAME\n" >> $LOG_FILE # for qsub
echo -e "\nScript file contents:\n" >> $LOG_FILE
cat "$0" >> $LOG_FILE
# ~~~~~~~~~~~~ #


# ~~~~~~ run FastQC ~~~~~~ #
# get the file name
# z=$(basename "$FASTQ1" )

# strip the extension
# z=${z/.fastq.gz*/} 

# set the outdir
# Outdir="$FASTQC_DIR/$z"
# mkdir -p "$Outdir"

fastqc --threads "$THREADS" --nogroup --outdir "$FASTQC_DIR" "$FASTQ1"
