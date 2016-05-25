#!/bin/bash
shopt -s extglob  #Enables extglob

# this script copies the results of a chip-seq analysis over to a Results dir for a client. 
# use different copy methods for different sets of data

# call script like this:
# scripts-copy-chipseq.sh {input/project dir full paths} {results dir full path}
# EXMPLE:
# ~/code/scripts-copy-chipseq_custom1.sh ~/projects/ ~/chip-seq/


#
##
### COPY THE BIG WIGS
if [ -d ${1}*/peaks/results ]; then
# make dir
mkdir -p $2alignment/bigwigs

# find the bigwig files
BIG_WIGS=$(ls ${1}*/alignments/results/*/*.bw)
for i in $BIG_WIGS; do
  # get the sample name from its dir name
  q=$(basename $(dirname $i ) );
  
  # strip out things we don't need in the name
  z=${q/*./};
  
  # copy the file with new name (don't overwrite if already present!)
  cp -an $i $2alignment/bigwigs/$z.bw
done
fi
###
##
#



#
##
### COPY THE PEAKS
if [ -d ${1}*/peaks/results ]; then
# make the dirs
mkdir -p $2peaks/{narrow,broad}

# all the files from the peak results that we will send to the client
ALL_PEAKS_FILES=$(ls ${1}*/peaks/results/peaks.macs_*/*/!(*@(.err|.id|.out|.sh)))
for i in $ALL_PEAKS_FILES; do
  # set broad/narrow based on the peak dir
  if [[ $i == *macs_narrow* ]]; then
    PEAKS_nar_brd="narrow"
  else
  
    if [[ $i == *macs_broad* ]]; then
      PEAKS_nar_brd="broad"
    else
      :
    fi
  fi
  
  # get the base dir name
  q=$(basename $(dirname $i ) );
  
  # fix the dir name
  z=${q/*./};
  
  # use the name as part of the filename
  beautiful_file_name=${z}_$(basename $i)

  # copy the file to the proper dir with the new name
  cp -an $i $2/peaks/$PEAKS_nar_brd/$beautiful_file_name

done
fi
###
##
#


#
##
### copy the PCA
if [ -d ${1}*/pca/results ]; then
  rsync -a --exclude "logs" --exclude 'job*' --exclude '*win*' ${1}*/pca/results/ $2pca
fi
# copy the Peaks table
if [ -d ${1}*/peaktable/results ]; then
  rsync -a --exclude "job*" --exclude "ref*" ${1}*/peaktable/results/peaktable.standard/ $2peaktable
fi
###
##
#

### PERMISSIONS
# allow group users to read and write
chmod g+rw -R $2*
chgrp -Rf results $2* # doesn't work yet..? Fail silently


shopt -u extglob # turn off the extglob
