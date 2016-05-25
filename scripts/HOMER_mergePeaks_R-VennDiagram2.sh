#!/bin/bash
# BED files with the peaks to overlap
tmp_outH3K4ME3="peaks_H3K4ME3.bed"
tmp_outH3K27AC="peaks_H3K27AC.bed"
# a sample ID
tmp_sampleID="ABC"

# HOMER mergePeaks
mergePeaks "$tmp_outH3K4ME3" "$tmp_outH3K27AC" -prefix mergepeaks -venn mergepeaks_venn

# the mergePeaks file outputs names:
tmp_mergeH3K4ME3="mergepeaks_${tmp_outH3K4ME3}"
tmp_mergeH3K27AC="mergepeaks_${tmp_outH3K27AC}"

# count the number of unique peaks
num_H3K4ME3=$(tail -n +2 $tmp_mergeH3K4ME3 | wc -l)
echo "num_H3K4ME3 is $num_H3K4ME3"
num_H3K27AC=$(tail -n +2 $tmp_mergeH3K27AC | wc -l)
echo "num_H3K27AC is $num_H3K27AC"

# count the number of peaks in common
num_overlap=$(tail -n +2 "mergepeaks_${tmp_outH3K4ME3}_${tmp_outH3K27AC}" | wc -l)

# plot the values in a pairwise venn in R
# # make sure the correct version of R is loaded:
module unload r
module load r/3.2.0
Rscript --slave --no-save --no-restore - "$tmp_sampleID" "$num_H3K4ME3" "$num_H3K27AC" "$num_overlap" <<EOF
  ## R code
  # load packages
  library('VennDiagram')
  library('gridExtra')
  # get script args, print them to console
  args <- commandArgs(TRUE); cat("Script args are:\n"); args
  SampleID<-args[1]
  peaks_H3K4ME3<-as.numeric(args[2])
  peaks_H3K27AC<-as.numeric(args[3])
  peaks_overlap<-as.numeric(args[4])
  # get filename for the plot PDF
  plot_filename<-paste0(SampleID,"_peaks.pdf") 
  # make a Venn object, don't print it yet
  venn<-draw.pairwise.venn(area1=peaks_H3K4ME3+peaks_overlap,area2=peaks_H3K27AC+peaks_overlap,cross.area=peaks_overlap,category=c('H3K4ME3','H3K27AC'),fill=c('red','blue'),alpha=c(0.3,0.3),cex=c(2,2,2),cat.cex=c(1.25,1.25),main=SampleID,ind=FALSE)
  # print it inside a PDF file, with a title
  pdf(plot_filename,width = 8,height = 8)
  grid.arrange(gTree(children=venn), top=SampleID) #, bottom="subtitle")
  dev.off()
EOF
