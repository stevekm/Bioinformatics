#!/bin/bash
# called by the qsub_motif_analysis.sh script
# USAGE: alignment_summary_stats_plots.sh /path/to/project/chipseq-analysis /path/to/project/chipseq-analysis_output genome
# $0 path to script
# $1 /path/to/project/chipseq-analysis
# $2 /path/to/output_dir 
# $3 <mm9/mm10/hg19/hg18> # pick one

module load homer/v4.6
# 
# echo "Inputs from here $1"
# echo "Outputs to here $2"

# get the number of threads from the qsub argument
THREADS=$NSLOTS 

# project dir
PROJ_INPUT_DIR=${1}

# output dir
MOTIF_OUT_DIR=${2}/motif_analysis

# get the refgenome
GENOME=$3

# set a dir for the refgenome parsing
PREPARSED_DIR=~/software/homer/${GENOME}/preparsed
mkdir -p $PREPARSED_DIR

# get the n'th sample for motif analysis
SAMPLES_INPUT=$(find ${PROJ_INPUT_DIR}/peaks/SOX9* -type f -name "summits.bed" | tr '\n' , | cut -d , -f $SGE_TASK_ID )

echo $SAMPLES_INPUT

# get the sample name from its dir
SAMPLE_NAME=$( dirname ${SAMPLES_INPUT} )
SAMPLE_NAME=$( basename ${SAMPLE_NAME} )
echo $SAMPLE_NAME

# make output dir for each sample
SAMPLE_OUT_DIR=$MOTIF_OUT_DIR/$SAMPLE_NAME
mkdir -p $SAMPLE_OUT_DIR

# run the motif analysis
echo "findMotifsGenome.pl $SAMPLES_INPUT $GENOME $SAMPLE_OUT_DIR -size 200 -preparsedDir $PREPARSED_DIR -p $THREADS"
findMotifsGenome.pl $SAMPLES_INPUT $GENOME $SAMPLE_OUT_DIR -size 200 -preparsedDir $PREPARSED_DIR -p $THREADS
