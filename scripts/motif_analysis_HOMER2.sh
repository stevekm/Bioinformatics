#!/bin/bash

##
## USAGE: motif_analysis_HOMER2.sh /path/to/project/motif-analysis_output /path/to/project/chipseq/sample.bed /path/to/params_file
##

# $0 path to script
# $1 /path/to/project/chipseq-analysis_dir
# $2 /path/to/output_dir 
# $3 /path/to/params_file # <mm9/mm10/hg19/hg18> # pick one

# set up logfile
LOG_FILE=$(basename $0).$(date -u +%Y%m%dt%H%M%S).${HOSTNAME:-$(hostname)}.$$.${RANDOM}.log

# process command-line inputs
# # if not enough args, output USAGE and exit
if (($# != 3)); then
  grep '^##' $0
  exit
fi


# # Input Script Args
MOTIF_OUT_DIR=$1
SAMPLE_INPUT_BED=$2
PARAMS=$3 

# # Get the params
source $PARAMS

# # Setup
# get the number of threads 
if [ $NSLOTS ]; then
  THREADS=$NSLOTS
  echo -e "NSLOTS is set to $NSLOTS, use for THREADS"
fi

# make sure the refgenome preparsed dir exists (from params)
mkdir -p $PREPARSED_DIR

# get the sample name from its dir
SAMPLE_NAME=$( dirname ${SAMPLE_INPUT_BED} )
SAMPLE_NAME=$( basename ${SAMPLE_NAME} )
echo -e "\nSample name is $SAMPLE_NAME"

# make output dir for each sample
SAMPLE_OUT_DIR=$MOTIF_OUT_DIR/$SAMPLE_NAME
mkdir -p $SAMPLE_OUT_DIR
cd $SAMPLE_OUT_DIR # HOMER makes tmp files so switch to here

# run the motif analysis on the sample file
# # only run the motif analysis if it hasn't already been run;
# # test if the dir has files in it already
if [[ -n $(find $SAMPLE_OUT_DIR -type f) ]]
then
  echo -e "\n OUTPUT DIR NOT EMPTY, SKIPPING MOTIF ANALYSIS \n"
else
  echo -e "\nRunning Motif Analysis\n"
  echo -e "Input file is $SAMPLE_INPUT_BED"
  echo -e "Genome is $GENOME"
  echo -e "Output dir is $SAMPLE_OUT_DIR"
  echo -e "Command is: \n"
  echo -e "findMotifsGenome.pl $SAMPLE_INPUT_BED $GENOME $SAMPLE_OUT_DIR -size $REGION_SIZE -preparsedDir $PREPARSED_DIR -p $THREADS \n"
  findMotifsGenome.pl $SAMPLE_INPUT_BED $GENOME $SAMPLE_OUT_DIR -size $REGION_SIZE -preparsedDir $PREPARSED_DIR -p $THREADS
  echo -e "findMotifsGenome.pl $SAMPLE_INPUT_BED $GENOME $SAMPLE_OUT_DIR -size $REGION_SIZE -preparsedDir $PREPARSED_DIR -p $THREADS" >> $SAMPLE_OUT_DIR/$LOG_FILE
fi

# Finding Instances and Locations of Specific Motifs
# # http://homer.salk.edu/homer/microarray/index.html
# # http://homer.salk.edu/homer/ngs/quantification.html
# # To recover the motif locations, you must first select the motifs you're interested in by getting the "motif file" output by HOMER.  
# # You can combine multiple motifs in single file if you like to form a "motif library".

# # create Motif Library from the known .motif files
MOTIF_LIB_FILE=$SAMPLE_OUT_DIR/known_motif_library.txt
rm -f $MOTIF_LIB_FILE # delte if already exists..
cat $SAMPLE_OUT_DIR/knownResults/*.motif > $MOTIF_LIB_FILE
# # I wrote this to add newlines between each motif but HOMER doesnt like them...
# rm -f $MOTIF_LIB_FILE
# for i in $(find $SAMPLE_OUT_DIR/knownResults -type f -name "*.motif"); do 
#   cat $i >> $MOTIF_LIB_FILE
#   echo "" >> $MOTIF_LIB_FILE
# done

# # find the motifs from the library
# # # check if the file already exists first
MOTIF_LOC_FILE=$SAMPLE_OUT_DIR/known_motif_locations.txt
rm -f $MOTIF_LOC_FILE # for debugging
if [[ -f $MOTIF_LOC_FILE ]]
then
  echo -e "MOTIF LOCATIONS ALREADY FOUND, SKIPPING STEP"
else
  echo -e "Finding motif locations"
  echo -e "\tInput: $MOTIF_LIB_FILE"
  echo -e "\tOutput: $MOTIF_LOC_FILE"
  annotatePeaks.pl tss $GENOME -size -300,50 -annStats $SAMPLE_OUT_DIR/annotation_stats.txt -m $MOTIF_LIB_FILE > $MOTIF_LOC_FILE
#   for i in $(find $SAMPLE_OUT_DIR/knownResults -type f -name "*.motif"); do 
#     annotatePeaks.pl tss $GENOME -size -300,50 -annStats $SAMPLE_OUT_DIR/annotation_stats_$(basename $i).txt -m $i >> $MOTIF_LOC_FILE
#   done
fi

# run the R script for parsing the motif table; set in params
$MOTIF_LOC_TABLE_SCRIPT $MOTIF_LOC_FILE

# ~~~~~~ save this as the separate params file; motif_analysis_HOMER2_params.sh
#!/bin/bash

# module load homer/v4.6

# set a dir for the refgenome parsing
# GENOME=mm9 # should contain the genome; mm9, mm10, etc.
# to find out which genomes are installed, check 
# ls /local/apps/homer/v4.6/data/genomes/

# set a dir to hold HOMER parsed refgenome data
# PREPARSED_DIR=/ifs/home/$(whoami)/software/homer/preparsed/${GENOME}
# # by defualt keep this hardlinked dir for preparesed

# Set the number of threads to be used; overridden by $NSLOTS if qsub is used
# THREADS=8

# location of the motif location table script
# MOTIF_LOC_TABLE_SCRIPT=/ifs/home/$(whoami)/projects/code/motif_analysis_HOMER2_LocTable.R

#  Selecting the size of the region for motif finding (-size # or -size given, default: 200)
# REGION_SIZE=200


# ~~~~~~ test things here

# project dir
# PROJ_INPUT_DIR=~/projects/
# 
# # output dir
# MOTIF_OUT_DIR=~/projects/motif_analysis
# 
# # params 
# PARAMS=~/projects/code/motif_analysis_HOMER2_params.sh
# 
# MOTIF_SCRIPT=~/projects/code/motif_analysis_HOMER2.sh
# 
# # Input_Entries=$(find ${PROJ_INPUT_DIR}/peaks/results/by-sample/peaks.macs_narrow/Nkx2* -type f -name "macs_summits.bed" | wc -l)
# 
# for i in $(find ${PROJ_INPUT_DIR}/peaks/results/by-sample/peaks.macs_narrow/Nkx2* -type f -name "macs_summits.bed"); do
# qsub -pe threaded 4-12 $MOTIF_SCRIPT $MOTIF_OUT_DIR $i $PARAMS;
# done

# 
# ~~~~~~ test things here
