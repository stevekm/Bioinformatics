#!/bin/bash

##
## USAGE: scripts-fastqc.sh /path/to/outdir /path/to/input_file
## This script will run FastQC; this script needs to be submitted with qsub
## 

# make sure that the correct number of script arguments were used
if [ $# != 2 ] # if not enough args provided
then
  grep '^##' $0 # print out lines from the script that start with '##' and exit
  exit
fi



# get the arguments
FASTQC_DIR="$1" # outdir
FASTQ1="$2" # input file
THREADS=$NSLOTS #  from qsub

# make the outdir
q=$(basename $(dirname "$FASTQ1" ) ) # the sample name
z=$(basename "$FASTQ1" ) # the file name
z=${z/.fastq.gz*/} # remove everything after 
Outdir="$FASTQC_DIR/$q/$z"
# echo  "$FASTQC_DIR/$z"
mkdir -p "$Outdir"

module load fastqc/0.11.4
fastqc --threads "$THREADS" --nogroup --outdir "$Outdir" "$FASTQ1"
# echo fastqc --threads "$THREADS" --nogroup --outdir "$FASTQC_DIR/$z" "$FASTQ1"
