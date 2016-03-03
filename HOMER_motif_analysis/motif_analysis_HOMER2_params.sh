#!/bin/bash

module load homer/v4.6

# set a dir for the refgenome parsing
GENOME=mm9 # should contain the genome; mm9, mm10, etc.
# to find out which genomes are installed, check 
# ls /local/apps/homer/v4.6/data/genomes/

# set a dir to hold HOMER parsed refgenome data
PREPARSED_DIR=/ifs/home/$(whoami)/software/homer/preparsed/${GENOME}
# # by defualt keep this hardlinked dir for preparesed

# Set the number of threads to be used; overridden by $NSLOTS if qsub is used
THREADS=8

# location of the motif location table script # this file will be updated and included later
MOTIF_LOC_TABLE_SCRIPT=/ifs/home/$(whoami)/projects/code/motif_analysis_HOMER2_LocTable.R
Peak_Types_RScript="/ifs/home/$(whoami)/projects/code/peak_types_plot.R"

#  Selecting the size of the region for motif finding (-size # or -size given, default: 200)
REGION_SIZE=200
