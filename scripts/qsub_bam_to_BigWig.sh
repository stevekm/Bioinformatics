#!/bin/bash

# this contains the qsub script that submits a job to the cluster

# directory containing fastq.gz files for paired end data
Input_Dir=~/project/input/fastq_dir

# the number of samples to be processed
Input_Entries=$(ls $Input_Dir/*R1.fastq.gz | wc -l)

# genome version to use
GENOME=mm10

# project out dir containing the aligned bam's to turn into BIGWIGs
BAM_INPUT_DIR=~/project/output/bam_dir

# dir to store the BIGWIGs in
BIGWIG_OUTPUT_DIR=~/project/output/bam_to_bigWig_py

# submit a job to the clutser
qsub -t 1-$Input_Entries -pe threaded 1 bam_to_BigWig.sh $GENOME $BAM_INPUT_DIR $BIGWIG_OUTPUT_DIR

# array job; 1 - number of samples
# a range of threads can be used for multithreaded tasks but this runs fine
# with just one thread for each process
