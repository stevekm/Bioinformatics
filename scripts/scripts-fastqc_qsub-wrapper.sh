#!/bin/bash

##
## USAGE: scripts-fastqc_qsub-wrapper.sh /path/to/outdir /path/to/project_dir
## This script will run FastQC on pipeline input fastq.gz files; this script needs to be submitted with qsub
## 

# make sure that the correct number of script arguments were used
if [ $# != 2 ] # if not enough args provided
then
  grep '^##' $0 # print out lines from the script that start with '##' and exit
  exit
fi

#~~~ get the arguments ~~~~#
FastQCscript="$(dirname "$0")/scripts-fastqc_qsub.sh" # script to be submitted with qsub

OutDir="$(readlink -m "$1")" # allow outdir that doesn't exist yet!

ProjDir="$2"

FastqDir="$(readlink -m "$2"/pipeline/inputs/fastq)" # parent directory with subdirectories per sample, each containing fastq.gz files

#~~~ Shell Options ~~~~#
# remember whether extglob was originally set, so we know whether to unset it
shopt -q extglob; extglob_set=$?
# set extglob if it wasn't originally set.
((extglob_set)) && shopt -s extglob
# Note, 0 (true) from shopt -q is "false" in a math context.
shopt -q nullglob; nullglob_set=$?
((nullglob_set)) && shopt -s nullglob
shopt -q globstar; globstar_set=$?
((globstar_set)) && shopt -s globstar
#~~~~~~~~~~~~~~~~~~#


# get the files to be procesed
FastqFiles=($FastqDir/**/*.fastq.gz)

for i in "${FastqFiles[@]}"; do
  echo "$i"
  qsub -pe threaded 1-8 "$FastQCscript" "$OutDir" "$i"
done

#~~~~~~~~~~~~~~~~~~#
# unset globs if it wasn't originally set
((extglob_set)) && shopt -u extglob
((nullglob_set)) && shopt -u nullglob
((globstar_set)) && shopt -u globstar
#~~~~~~~~~~~~~~~~~~#
# how many files are there??
# ArrayJobLength=$(echo "${FastqFiles[@]}" | tr " " "\n" | wc -l)
# don't use this right now; trouble passing array object to array job.. 
