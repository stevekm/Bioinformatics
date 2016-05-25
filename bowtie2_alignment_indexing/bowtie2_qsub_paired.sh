#!/bin/bash
# print the script commands to stderr on execution; stderr gets saved in the qsub log files
set -x
##
## USAGE: bowtie2_qsub_paired.sh /path/to/outdir /path/to/input_file_R1.fastq.gz /path/to/input_file_R2.fastq.gz /path/to/ref_genome 
## This script will run bowtie2; this script should be submitted with qsub!
## This version takes paired end reads

# ~~~~~~ permissions ~~~~~~ #
# To make stuff group-writeable
umask 007

# ~~~~~~ script args processing ~~~~~~ #
# # if not enough args, output USAGE info lines (which start with '##') and exit
if (($# != 4)); then
  grep '^##' $0
  exit
fi

# ~~~~~~ bowtie params ~~~~~~ #
# these need to be installed on the cluster
module load bowtie2 # http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
module load samtools/1.2.1
mapq=30

# get the script arguments
OUTDIR="$1" # outdir
FASTQ_R1="$2" # Read 1 input file
FASTQ_R2="$3" # Read 2 file
GENOME="$4" # the name of the genome to be used e.g. hg19, mm10, etc
# THREADS is defined by qsub; if not, set with 8 threads
THREADS=${NSLOTS:=8} 
echo "OUTDIR is $OUTDIR"
echo "FASTQ_R1 is $FASTQ_R1"
echo "FASTQ_R2 is $FASTQ_R2"
echo "GENOME is $GENOME"
echo "THREADS is $THREADS"

# ~~~~~~ Script Logging ~~~~~~~ #
# save a hardcopy of this script file in case there are future changes
zz=$(basename $0)
# get the script dir
za=$(dirname $0)
# for regular script running (no qsub)
# LOG_FILE=log.$(basename $0).$(date -u +%Y%m%dt%H%M%S).${HOSTNAME:-$(hostname)}.$$.${RANDOM}
# for qsub script submission
LOG_FILE=scriptlog.$JOB_NAME.$(date -u +%Y%m%dt%H%M%S).${HOSTNAME:-$(hostname)}.$$.${RANDOM}
echo "This is the log file for the script." > $LOG_FILE
echo -e "\n pwd is:\n$PWD\n" >> $LOG_FILE
# echo -e "\nScript file is:\n$0\n" >> $LOG_FILE # for regular script usage (no qsub)
echo -e "\nScript file is:\n$JOB_NAME\n" >> $LOG_FILE # for qsub
echo -e "\nScript file contents:\n" >> $LOG_FILE
cat "$0" >> $LOG_FILE
# ~~~~~~~~~~~~ #

# ~~~~~~ run command ~~~~~~ #
# align with bowtie2
# bowtie aligns, output sam default
# pipe to samtools view, optional mapping quality score filter mapq, convert to bam
# pipe to samtools sort to sort the entries

tmpOUTFILE="$(basename "$FASTQ_R1")"
tmpOUTFILE="${tmpOUTFILE//.fastq.gz/}"

# location of the bowtie2 genome
# /local/data/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome
tmpGenome="/local/data/iGenomes/Homo_sapiens/UCSC/${GENOME}/Sequence/Bowtie2Index/genome"
echo "bowtie2 --threads $THREADS --local -x $tmpGenome -q -1 $FASTQ_R1 -2 $FASTQ_R2 | samtools view -@ $THREADS -Sb1 - | samtools sort -m 10G -@ $THREADS - $tmpOUTFILE"
bowtie2 --threads "$THREADS" --local -x "$tmpGenome" -q -1 "$FASTQ_R1" -2 "$FASTQ_R2" | samtools view -@ "$THREADS" -Sb1 - | samtools sort -m 10G -@ "$THREADS" - "$tmpOUTFILE"
samtools index "${tmpOUTFILE}.bam"

# some user manuals that are useful for these tools:
# http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
# http://www.htslib.org/doc/samtools-0.1.19.html
# http://samtools.github.io/hts-specs/SAMv1.pdf
# https://edwards.sdsu.edu/research/creating-indexed-bam-files-from-bowtie-alignments/

