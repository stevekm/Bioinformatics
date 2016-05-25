#!/bin/bash
set -x # debugging output to stderr

# load the correct version of R
module unload r; module load r/3.3.0

# make sure the script is executable
chmod +x multi_peaks_UpSet_plot.R

# run the script to make an UpSet plot in R based on the HOMER venn.txt
Rscript --vanilla multi_peaks_UpSet_plot.R "All Samples" venn.txt
