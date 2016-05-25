#!/bin/bash
# USAGE: HOMER_KnownMotif_gene_name_extraction.sh /path/to/motif_analysis_dir
# this script will parse all .motif files in the HOMER knownResults and homerResults dirs created by HOMER
# and get the names of the genes for which motifs were found
# output will be saved as a text list of the gene names
# # note: some items will not actually be gene names! Make sure to check this

# get the motif_analysis_dir
MOTIF_DIR=$1
echo $MOTIF_DIR

# output the Gene names to a new text file
cat $MOTIF_DIR/knownResults/*.motif | grep '>' | cut -f2 | cut -d '(' -f1 > $MOTIF_DIR/knownResults_gene_list.txt
