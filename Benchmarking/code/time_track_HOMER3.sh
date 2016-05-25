#!/bin/bash

##
## USAGE: time_track_HOMER3.sh /path/to/project/motif-analysis_output_dir /path/to/project/chipseq/sample.bed /path/to/timer_log_file
## This script will run HOMER motif analysis on a bed file, tracking time to completion

# ~~~~~~ script args ~~~~~~ #
# # if not enough args, output USAGE info lines (which start with '##') and exit
if (($# != 3)); then
  grep '^##' $0
  exit
fi

echo -e "start\t${NSLOTS}\t$(date +%s)" > start.txt
echo -e "start\t${NSLOTS}\t$(date +%s)"
echo "Log file is $TIMER_LOGFILE"

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


# # Input Script Args
OUT_DIR="$1" # should be the pwd as per qsub params
SAMPLE_INPUT_BED="$2" 
TIMER_LOGFILE="$3" # file to record the times in
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



for (( i=1; i<=32; i+=1 )); do
  echo "$i"
  THREADS="$i"
  tmp_Subdir="${OUT_DIR}/$i"
  mkdir -p "$tmp_Subdir"
  # cd "$tmp_Subdir"
  
  # run the motif analysis
  start_time1=`date +%s`
  echo -e "pwd is $(pwd)\n\nfindMotifsGenome.pl $SAMPLE_INPUT_BED $GENOME $OUT_DIR -size $REGION_SIZE -preparsedDir $PREPARSED_DIR -p $THREADS && echo -e $THREADS\t$(expr $(date +%s) - $start_time1) >> $TIMER_LOGFILE"
  findMotifsGenome.pl "$SAMPLE_INPUT_BED" "$GENOME" "$tmp_Subdir" -size "$REGION_SIZE" -preparsedDir "$PREPARSED_DIR" -p "$THREADS" && echo -e "$THREADS\t$(expr $(date +%s) - $start_time1)" >> "$TIMER_LOGFILE"

done


echo -e "stop\t${NSLOTS}\t$(date +%s)"
echo -e "stop\t${NSLOTS}\t$(date +%s)" > stop.txt


