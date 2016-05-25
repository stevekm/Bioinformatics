#!/bin/bash

##
## USAGE: make_tagDirectory_HOMER-bed.sh /path/to/project/output_dir /path/to/project/chipseq/sample.bed /path/to/params_file
## This script will run HOMER makeTagDirectory on a .BED file


# ~~~~~~ script args ~~~~~~ #

# $0 path to script, except in qsub submission.. 
# $1 /path/to/project/motif-analysis_output
# $2 /path/to/project/chipseq/sample.bed
# $3 /path/to/params_file # <mm9/mm10/hg19/hg18> # pick one

# # if not enough args, output USAGE info lines (which start with '##') and exit
if (($# != 3)); then
  grep '^##' $0
  exit
fi


# # Input Script Args
TAG_OUT_DIR="$1"
SAMPLE_INPUT_BED="$2"
PARAMS="$3"
source "$PARAMS"
# ~~~~~~~~~~~~ #


# ~~~~~~ Script Logging ~~~~~~~ #
# get the script file 
zz=$(basename $0)
# get the script dir
za=$(dirname $0)
# echo $0
# echo $zz
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

echo -e "\nParams file is:\n" >> $LOG_FILE
echo -e "$PARAMS" >> $LOG_FILE
echo -e "\nParams file contents:\n" >> $LOG_FILE
cat "$PARAMS" >> $LOG_FILE
# ~~~~~~~~~~~~ #


# # Setup
# get the number of threads 
if [ $NSLOTS ]; then
  THREADS=$NSLOTS
  echo -e "NSLOTS is set to $NSLOTS, use for THREADS"
else
  THREADS="8"
fi

# make sure the refgenome preparsed dir exists (from params)
mkdir -p $PREPARSED_DIR

# make the Tag Directory
echo -e "Command is:\n"
echo -e "makeTagDirectory $TAG_OUT_DIR -genome $GENOME -format bed $SAMPLE_INPUT_BED"
makeTagDirectory "$TAG_OUT_DIR" -genome "$GENOME" -format bed "$SAMPLE_INPUT_BED"
