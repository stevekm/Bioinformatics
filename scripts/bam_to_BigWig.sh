#!/bin/bash

# USAGE: bam_to_BigWig.sh $GENOME $BAM_INPUT_DIR $BIGWIG_OUTPUT_DIR
# use with qsub array job based on number of samples!! Like this:
########
# EXAMPLE
# the fastq.gz files are in this dir:
# Input_Dir=~/projects/RNAseq-2015-09-29/inputs 
# # count the number of entries to be processed:
# Input_Entries=$(ls $Input_Dir/*R1.fastq.gz | wc -l)
# 
# # the name of the referrence genome
# GENOME=mm9
# # the location of the bams to be turned into BigWigs:
# BAM_INPUT_DIR=/ifs/home/username/projects/RNAseq-2015-09-29/output_mm9/bam
# # the location to output the BigWigs to:
# BIGWIG_OUTPUT_DIR=/ifs/home/username/projects/RNAseq-2015-09-29/output_mm9/bam_to_bigWig_py
# # ~/projects/RNAseq-2015-09-29/code/bam_to_BigWig.sh $GENOME $BAM_INPUT_DIR $BIGWIG_OUTPUT_FILE
# qsub -t 1-$Input_Entries -pe threaded 1 ~/projects/RNAseq-2015-09-29/code/bam_to_BigWig.sh $GENOME $BAM_INPUT_DIR $BIGWIG_OUTPUT_DIR

# need these for the python functions to work

module load samtools/0.1.19
module load bedtools/2.22.0

# these are the arguments:
# echo -e "script dir is $0 \n"
# 
# echo -e "genome is $1 \n"
# 
# echo -e "bam input dir is $2 \n"
# 
# echo -e "bigwig output file is $3 \n"
# 
# echo "Task ID is ${SGE_TASK_ID}"
GENOME=$1
BAM_INPUT_DIR=$2
BIGWIG_OUTPUT_DIR=$3

# make sure the output dir exists
mkdir -p $BIGWIG_OUTPUT_DIR

# select the BAM_INPUT_FILE from the dir
BAM_INPUT_FILE=$( basename $(ls $BAM_INPUT_DIR/*.bam | tr '\n' , | cut -d , -f $SGE_TASK_ID) )
BAM_INPUT_FILE=${BAM_INPUT_DIR}/${BAM_INPUT_FILE}

echo -e "BAM Input File is $BAM_INPUT_FILE \n"

# set the output bigwig file name... 
BIGWIG_OUTPUT_FILE=$( basename $(ls $BAM_INPUT_DIR/*.bam | tr '\n' , | cut -d , -f $SGE_TASK_ID) )
BIGWIG_OUTPUT_FILE=$BIGWIG_OUTPUT_DIR/${BIGWIG_OUTPUT_FILE//.bam/}.bw
echo -e "BigWig Output File is $BIGWIG_OUTPUT_FILE \n"

echo "Calling Python..."
python bam_to_BigWig.py $GENOME $BAM_INPUT_FILE $BIGWIG_OUTPUT_FILE
