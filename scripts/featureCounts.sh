#!/bin/bash

# this script will run featureCounts with some pre-selected parameters for paired end data

# FROM THE MANUAL:
# Summarize multiple paired-end datasets:
# featureCounts -p -t exon -g gene_id -a annotation.gtf -o counts.txt library1.bam library2.bam library3.bam


GTF=$(echo /local/data/iGenomes/Mus_musculus/UCSC/${1}/Annotation/Genes/genes.gtf)
# THREADS=4
THREADS=$NSLOTS
# SGE_TASK_ID=1

# BAM_DIR=~/projects/output/bam
BAM_DIR=$2

# OUT_DIR=~/projects/output/feature_counts
OUT_DIR=$3

mkdir -p $OUT_DIR

# Get the name of the sample from the bam file
# SGE_TASK_ID is supplied by qsub array job ID during the script's qsub call
# # here, it tells the script which entry in the list of files to select
SAMPLE_NAME=$( basename $(ls $BAM_DIR/*.bam | tr '\n' , | cut -d , -f $SGE_TASK_ID) )
SAMPLE_NAME=${SAMPLE_NAME//.bam/}

##########

# featureCounts generates TONS of temp files in the current directory; 
# cd to the output dir and do the work there
cd $OUT_DIR

module load subread/1.4.6-p3
echo "featureCounts -p -T $THREADS -g gene_name -a $GTF -o $OUT_DIR/${SAMPLE_NAME}_feature_counts.txt $BAM_DIR/${SAMPLE_NAME}.bam"
featureCounts -p -T $THREADS -g gene_name -a $GTF -o $OUT_DIR/${SAMPLE_NAME}_feature_counts.txt $BAM_DIR/${SAMPLE_NAME}.bam

