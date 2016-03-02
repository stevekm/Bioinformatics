#!/bin/bash

##
## USAGE: motif_analysis_HOMER3.sh /path/to/project/motif-analysis_output /path/to/project/chipseq/sample.bed /path/to/params_file
## This script will run HOMER annotatePeaks and findMotifsGenome on a .BED file


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
MOTIF_OUT_DIR="$1"
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

echo -e "\nRunning Motif Analysis\n"
echo -e "Input file is $SAMPLE_INPUT_BED"
echo -e "Genome is $GENOME" # from params
echo -e "Output dir is $MOTIF_OUT_DIR"
echo -e "Command is: \n"
echo -e "findMotifsGenome.pl $SAMPLE_INPUT_BED $GENOME $MOTIF_OUT_DIR -size $REGION_SIZE -preparsedDir $PREPARSED_DIR -p $THREADS \n"
findMotifsGenome.pl $SAMPLE_INPUT_BED $GENOME $MOTIF_OUT_DIR -size $REGION_SIZE -preparsedDir $PREPARSED_DIR -p $THREADS

# Finding Instances and Locations of Specific Motifs
# # http://homer.salk.edu/homer/microarray/index.html
# # http://homer.salk.edu/homer/ngs/quantification.html
# # To recover the motif locations, you must first select the motifs you're interested in by getting the "motif file" output by HOMER.  
# # You can combine multiple motifs in single file if you like to form a "motif library".

# # create Motif Library from the known .motif files
MOTIF_LIB_FILE=$MOTIF_OUT_DIR/known_motif_library.txt
rm -f $MOTIF_LIB_FILE # delete if already exists.. for debugging purposes
# print all the motifs into one file
cat $MOTIF_OUT_DIR/knownResults/*.motif > $MOTIF_LIB_FILE

# # find the motifs from the library
# # # check if the file already exists first # consider removing this check & overwriting old outputs on re-runs of script
MOTIF_LOC_FILE=$MOTIF_OUT_DIR/known_motif_locations.txt
rm -f $MOTIF_LOC_FILE # for debugging
if [[ -f $MOTIF_LOC_FILE ]]
then
  echo -e "MOTIF LOCATIONS ALREADY FOUND, SKIPPING STEP"
else
  echo -e "Finding motif locations"
  echo -e "\tInput: $MOTIF_LIB_FILE"
  echo -e "\tOutput: $MOTIF_LOC_FILE"
  echo -e "Command is: \n"
  echo -e "annotatePeaks.pl tss $GENOME -size -300,50 -annStats $MOTIF_OUT_DIR/annotation_stats.txt -m $MOTIF_LIB_FILE > $MOTIF_LOC_FILE \n"
  annotatePeaks.pl tss $GENOME -size -300,50 -annStats $MOTIF_OUT_DIR/annotation_stats.txt -m $MOTIF_LIB_FILE > $MOTIF_LOC_FILE
fi

# run the R script for parsing the motif table; set in params
# $MOTIF_LOC_TABLE_SCRIPT $MOTIF_LOC_FILE
# TODO : work on this later; script is already in the parent repo but needs updating, etc.
