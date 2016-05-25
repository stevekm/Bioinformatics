#!/bin/bash
echo -e "start\t${NSLOTS}\t$(date +%s)" > start.txt
echo -e "start\t${NSLOTS}\t$(date +%s)"
##
## USAGE: time_track_HOMER2.sh /path/to/project/motif-analysis_output /path/to/project/chipseq/sample.bed /path/to/timer_log_file
## This script will run HOMER motif analysis on a bed file, tracking time to completion

# ~~~~~~ HOMER Params ~~~~~~ #
module load homer/v4.6
module unload r
module load r/3.2.0

# set a dir for the refgenome parsing
GENOME=mm9 # should contain the genome; mm9, mm10, etc.

# set a dir to hold HOMER parsed refgenome data
PREPARSED_DIR=/ifs/home/$(whoami)/software/homer/preparsed/${GENOME}
mkdir -p $PREPARSED_DIR

#  Selecting the size of the region for motif finding (-size # or -size given, default: 200)
REGION_SIZE=200
# ~~~~~~ # ~~~~~~ # ~~~~~~ #


# ~~~~~~ script args ~~~~~~ #

# $0 path to script, except in qsub submission.. 
# $1 /path/to/project/motif-analysis_output
# $2 /path/to/project/chipseq/sample.bed

# # if not enough args, output USAGE info lines (which start with '##') and exit
if (($# != 3)); then
  grep '^##' $0
  exit
fi


# # Input Script Args
MOTIF_OUT_DIR="$1" # should be the pwd as per qsub params
SAMPLE_INPUT_BED="$2"
TIMER_LOGFILE="$3"
THREADS="$NSLOTS"
# ~~~~~~~~~~~~ #


# ~~~~~~ Script Logging ~~~~~~~ #
# get the script file 
zz=$(basename $0)
# get the script dir
za=$(dirname $0)
# for regular script running (no qsub)
# LOG_FILE=log.$(basename $0).$(date -u +%Y%m%dt%H%M%S).${HOSTNAME:-$(hostname)}.$$.${RANDOM}
# for qsub script submission
# LOG_FILE=log.$JOB_NAME.$(date -u +%Y%m%dt%H%M%S).${HOSTNAME:-$(hostname)}.$$.${RANDOM}
LOG_FILE=log.$JOB_NAME.${Script_start_time1}.${HOSTNAME:-$(hostname)}.$$.${RANDOM}
echo "This is the log file for the script." > $LOG_FILE
echo -e "\n pwd is:\n$PWD\n" >> $LOG_FILE
# echo -e "\nScript file is:\n$0\n" >> $LOG_FILE # for regular script usage (no qsub)
echo -e "\nScript file is:\n$JOB_NAME\n" >> $LOG_FILE # for qsub
echo -e "\nScript file contents:\n" >> $LOG_FILE
cat "$0" >> $LOG_FILE
# ~~~~~~~~~~~~ #

# run the motif analysis
start_time1=`date +%s`
echo -e "pwd is $(pwd)\n\nfindMotifsGenome.pl \n$SAMPLE_INPUT_BED \n$GENOME \n$MOTIF_OUT_DIR \n-size $REGION_SIZE \n-preparsedDir $PREPARSED_DIR \n-p $THREADS \n&& \necho \n-e $THREADS\t$(expr $(date +%s) - $start_time1) \n>> \n$TIMER_LOGFILE"

findMotifsGenome.pl "$SAMPLE_INPUT_BED" "$GENOME" "$MOTIF_OUT_DIR" -size "$REGION_SIZE" -preparsedDir "$PREPARSED_DIR" -p "$THREADS" && echo -e "$THREADS\t$(expr $(date +%s) - $start_time1)\t$JOB_ID" >> "$TIMER_LOGFILE"




echo -e "stop\t${NSLOTS}\t$(date +%s)"
echo -e "stop\t${NSLOTS}\t$(date +%s)" > stop.txt


