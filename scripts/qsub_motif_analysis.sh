#!/bin/bash
## USAGE qsub_motif_analysis.sh /path/to/project/chipseq-analysis /path/to/project/chipseq-analysis_output genome


# project dir
PROJ_INPUT_DIR=${1}

# output dir
MOTIF_OUT_DIR=${2}

# get the refgenome
GENOME=${3}

# the location of the motif analysis script
MOTIF_SCRIPT=/path/to/code/motif_analysis_HOMER.sh

# get the number of samples to be processed
Input_Entries=$(find ${PROJ_INPUT_DIR}/peaks/SOX9* -type f -name "summits.bed" | wc -l)

qsub -t 1-$Input_Entries -pe threaded 4-32 $MOTIF_SCRIPT $PROJ_INPUT_DIR $MOTIF_OUT_DIR $GENOME
