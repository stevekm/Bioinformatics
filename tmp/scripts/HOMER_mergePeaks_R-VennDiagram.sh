#!/bin/bash

# this is a workflow script to use as template for setting up HOMER mergePeaks and R VennDiagram

# for both MACS peaks by sample, group
for q in sample group; do
tmp_type="$q"

# load HOMER on the HPC session
module load homer/v4.6

# project directory
ProjDir="$HOME/projects/SmithLab_ChIpSeq_2016-03-10"

# the peaks files are in subdirectories burried in here:
peaks_resultsDir="${ProjDir}/pipeline/peaks/results"

# place to put the output
testDir="${ProjDir}/project_notes/peak_overlap/macs_broad_by_${tmp_type}"

# samplesheet with the data we need for sample ID's and such
samplesheet="${ProjDir}/inputs/sample-sheet.tsv"

# make the dir's and switch to them
mkdir -p "$testDir"
cd "$testDir"

# read in the samplesheet without header, cut for field 9 which has a list of the unique Sample ID's, iterate over them
for i in $(tail -n +2 $samplesheet | cut -f9 | sort -u); do
  # get the sample ID, setup the dir
  echo  "$i"
  tmp_sampleID="$i"
  
  tmp_outDir="${testDir}/${tmp_sampleID}"
  mkdir -p "$tmp_outDir"
  cd "$tmp_outDir"
  
  # find the peak files for each histone mark, copy to pwd because thats easier to work with
  tmp_bedfile_H3K27AC="${tmp_outDir}/${tmp_sampleID}_R_H3K27AC_peaks.bed"
  echo "tmp_bedfile_H3K27AC is $tmp_bedfile_H3K27AC"
  
  tmp_bedfile_H3K4ME3="${tmp_outDir}/${tmp_sampleID}_R_H3K4ME3_peaks.bed"
  echo "tmp_bedfile_H3K27AC is $tmp_bedfile_H3K4ME3"
  
  find "$peaks_resultsDir" -type f -name "peaks.bed" -path "*/peaks.by_${tmp_type}.macs_broad/*" -path "*/${tmp_sampleID}-R-H3K27AC*" -exec cp {} "${tmp_bedfile_H3K27AC}" \;
  
  find "$peaks_resultsDir" -type f -name "peaks.bed" -path "*/peaks.by_${tmp_type}.macs_broad/*" -path "*/${tmp_sampleID}-R-H3K4ME3*" -exec cp {} "${tmp_bedfile_H3K4ME3}" \;
  
  # make sure BOTH files exist before proceeding..
  if [[ -f ${tmp_bedfile_H3K4ME3} && -f ${tmp_bedfile_H3K27AC} ]]; then
    echo "Both files found"
    
    tmp_outH3K4ME3="$(basename ${tmp_bedfile_H3K4ME3})"
    tmp_outH3K27AC="$(basename ${tmp_bedfile_H3K27AC})"
    
    # merge the peak files with HOMER
    # do not pass file paths to mergePeaks! It craps out if you do, only pass pwd files
    mergePeaks "$tmp_outH3K4ME3" "$tmp_outH3K27AC" -prefix mergepeaks -venn mergepeaks_venn
    
    # the mergePeaks outputs names
    tmp_mergeH3K4ME3="mergepeaks_${tmp_outH3K4ME3}"
    tmp_mergeH3K27AC="mergepeaks_${tmp_outH3K27AC}"
      
    # get the number of peaks in each file; this is easier than parsing the venn.txt file...
    # unique peaks
    num_H3K4ME3=$(tail -n +2 $tmp_mergeH3K4ME3 | wc -l)
    echo "num_H3K4ME3 is $num_H3K4ME3"
    num_H3K27AC=$(tail -n +2 $tmp_mergeH3K27AC | wc -l)
    echo "num_H3K27AC is $num_H3K27AC"
    
    # peaks in common
    num_overlap=$(tail -n +2 "mergepeaks_${tmp_outH3K4ME3}_${tmp_outH3K27AC}" | wc -l)
    
    # load the correct R version on the HPC session
    module unload r
    module load r/3.2.0
    # need this version of R for VennDiagram pacakge
    # pass the commands to R in a heredoc from stdin (I love doing this)
    Rscript --slave --no-save --no-restore - "$tmp_sampleID" "$num_H3K4ME3" "$num_H3K27AC" "$num_overlap" <<EOF
      ## R code
      # load packages
      library('VennDiagram')
      library('gridExtra')
      # sanity check
      cat("\nR loaded\n")
      args <- commandArgs(TRUE); cat("Script args are:\n"); args
      # assign the args properly
      SampleID<-args[1]; peaks_H3K4ME3<-as.numeric(args[2]); peaks_H3K27AC<-as.numeric(args[3]); peaks_overlap<-as.numeric(args[4]); plot_filename<-paste0(SampleID,"_peaks.pdf") 
      # save the venn diagram object
      venn<-draw.pairwise.venn(area1=peaks_H3K4ME3+peaks_overlap,area2=peaks_H3K27AC+peaks_overlap,cross.area=peaks_overlap,category=c('H3K4ME3','H3K27AC'),fill=c('red','blue'),alpha=c(0.3,0.3),cex=c(2,2,2),cat.cex=c(1.25,1.25),main=SampleID,ind=FALSE)
      # open a PDF, and print the Venn with a title
      pdf(plot_filename,width = 8,height = 8)
      grid.arrange(gTree(children=venn), top=SampleID) #, bottom="subtitle")
      dev.off()
EOF
    else
    echo "A file is missing!"
  fi
  echo -e "\n\n"
  cd "$testDir"
done

# make a common dir for the PDFs and copy them there
mkdir -p "${testDir}/all_pdf"
find . -type f -name "*.pdf" -exec cp {} "${testDir}/all_pdf/" \;

done

