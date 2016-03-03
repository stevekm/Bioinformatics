#!/bin/bash

##
## USAGE: time_track_wrapper.sh /path/to/project/motif-analysis_output_dir /path/to/project/chipseq/sample.bed /path/to/script /path/to/timer_log_file
## This script will run HOMER motif analysis on a bed file repeatedly with a range of CPU threads from 32 to 1, tracking time to completion

# ~~~~~~ script args ~~~~~~ #
# # if not enough args, output USAGE info lines (which start with '##') and exit
if (($# != 4)); then
  grep '^##' $0
  exit
fi

# # Input Script Args
OUT_DIR="$1" # should be the pwd as per qsub params
SAMPLE_INPUT_BED="$2" 
tmpSCRIPT="$3" # script to get submitted
TIMER_LOGFILE="$4" # file to record the times in
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
LOG_FILE=log.$JOB_NAME.${Script_start_time1}.${HOSTNAME:-$(hostname)}.$$.${RANDOM} # writes in the current directory
echo "This is the log file for the script." > $LOG_FILE
echo -e "\n pwd is:\n$PWD\n" >> $LOG_FILE
# echo -e "\nScript file is:\n$0\n" >> $LOG_FILE # for regular script usage (no qsub)
echo -e "\nScript file is:\n$JOB_NAME\n" >> $LOG_FILE # for qsub
echo -e "\nScript file contents:\n" >> $LOG_FILE
cat "$0" >> $LOG_FILE
# ~~~~~~~~~~~~ #

# create the timer file
# echo -n "" > "$TIMER_LOGFILE"
echo "Log file is $TIMER_LOGFILE"

# set up the qsub submissions
# for i in {32..8}; do
for (( i=2; i<=32; i+=2 )); do
  echo "$i"
  THREADS="$i"
  echo -e "qsub -q all.q -wd ${OUT_DIR}/${i} \n-o :${OUT_DIR}/${i} \n-e :${OUT_DIR}/${i} \n-pe threaded $THREADS \n$tmpSCRIPT \n${OUT_DIR}/${i} \n${SAMPLE_INPUT_BED} \n${TIMER_LOGFILE}"
  qsub -q all.q -wd "${OUT_DIR}/${i}" -o :${OUT_DIR}/${i}/ -e :${OUT_DIR}/${i}/ -pe threaded "$THREADS" "$tmpSCRIPT" "${OUT_DIR}/${i}" "${SAMPLE_INPUT_BED}" "${TIMER_LOGFILE}"
done



 

