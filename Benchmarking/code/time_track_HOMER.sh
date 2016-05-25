#!/bin/bash
Script_start_time1=$(date -u +%Y%m%dt%H%M%S)
##
## USAGE: time_track_HOMER.sh /path/to/project/motif-analysis_output /path/to/project/chipseq/sample.bed
## This script will run HOMER motif analysis on a bed file repeatedly with a range of CPU threads from 32 to 1, tracking time to completion

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
if (($# != 2)); then
  grep '^##' $0
  exit
fi


# # Input Script Args
MOTIF_OUT_DIR="$1" # should be the pwd as per qsub params
SAMPLE_INPUT_BED="$2"

tmp_SAMPLE=$(basename $SAMPLE_INPUT_BED)
tmp_SAMPLE=${tmp_SAMPLE%.*}
TIMER_LOGFILE="${MOTIF_OUT_DIR}/timer_log_${tmp_SAMPLE}_${Script_start_time1}.txt"
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

# create the timer file
echo -n "" > "$TIMER_LOGFILE"
echo "Log file is $TIMER_LOGFILE"
# 
# # make a tmp outdir to wrangle the motif output easier # forget this, it gets messed up in qsub
# MOTIF_OUT_DIR_tmp="${MOTIF_OUT_DIR_tmp}/tmp"
# # delete it if it already exists
# if [ -d "$MOTIF_OUT_DIR_tmp" ]; then
#   rm -rf "$MOTIF_OUT_DIR_tmp"
# else
#   mkdir -p "$MOTIF_OUT_DIR_tmp"
# fi

for i in {32..1}; do
  THREADS="$i"
  start_time1=`date +%s`
  # echo -ne "$THREADS\t" >> $TIMER_LOGFILE && echo "$(expr $(date +%s) - $start_time1)" >> $TIMER_LOGFILE
  findMotifsGenome.pl "$SAMPLE_INPUT_BED" "$GENOME" "$MOTIF_OUT_DIR" -size "$REGION_SIZE" -preparsedDir "$PREPARSED_DIR" -p "$THREADS" && echo -e "$THREADS\t$(expr $(date +%s) - $start_time1)" >> "$TIMER_LOGFILE"
  
#   start_time2=`date +%s`
#   echo -ne "$THREADS\t" >> $TIMER_LOGFILE && echo "$(expr $(date +%s) - $start_time2)" >> $TIMER_LOGFILE
#   start_time3=`date +%s`
#   echo -ne "$THREADS\t" >> $TIMER_LOGFILE && echo "$(expr $(date +%s) - $start_time3)" >> $TIMER_LOGFILE
#   # <command-to-execute> && echo run time is $(expr `date +%s` - $start_time1) s
#   
  
  # echo -e "findMotifsGenome.pl $SAMPLE_INPUT_BED $GENOME $MOTIF_OUT_DIR -size $REGION_SIZE -preparsedDir $PREPARSED_DIR -p $THREADS \n"
  #findMotifsGenome.pl $SAMPLE_INPUT_BED $GENOME $MOTIF_OUT_DIR -size $REGION_SIZE -preparsedDir $PREPARSED_DIR -p $THREADS
done





 

