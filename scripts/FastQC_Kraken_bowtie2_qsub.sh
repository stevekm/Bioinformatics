#!/bin/bash
# print the script commands to stderr on execution; stderr gets saved in the qsub log files
set -x
##
## USAGE: FastQC_Kraken_bowtie2_qsub.sh /path/to/outdir /path/to/input_file_R1.fastq.gz /path/to/input_file_R2.fastq.gz /path/to/ref_genome 
## This script will run bowtie2; this script should be submitted with qsub!
## 

# ~~~~~~ permissions ~~~~~~ #
# To make stuff group-writeable
umask 007

# ~~~~~~ script args processing ~~~~~~ #
# # if not enough args, output USAGE info lines (which start with '##') and exit
if (($# != 5)); then
  grep '^##' $0
  exit
fi

# ~~~~~~ params ~~~~~~ #
# these need to be installed on the cluster
export KRAKEN_DEFAULT_DB="$HOME/ref/Kraken/nt-48G"
export KRAKEN_NUM_THREADS="$THREADS"
module load fastqc/0.11.4
module load bowtie2 # http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
module load samtools/1.2.1
mapq=30

# get the script arguments
OUTDIR="$1" # outdir
FASTQ_R1="$2" # Read 1 input file
FASTQ_R2="$3" # Read 2 file
GENOME="$4" # the name of the genome to be used e.g. hg19, mm10, etc
SAMPLEID="$5"
# THREADS is defined by qsub; if not, set with 8 threads
THREADS=${NSLOTS:=8} 
echo "OUTDIR is $OUTDIR"
echo "FASTQ_R1 is $FASTQ_R1"
echo "FASTQ_R2 is $FASTQ_R2"
echo "GENOME is $GENOME"
echo "THREADS is $THREADS"

# ~~~~~~ Script Logging ~~~~~~~ #
# save a hardcopy of this script file 
zz=$(basename $0)
# get the script dir
za=$(dirname $0)
# set the script log file, use JOBNAME if its set by qsub
LOG_FILE=scriptlog.${JOB_NAME:=$(basename $0)}.$(date -u +%Y%m%dt%H%M%S).${HOSTNAME:-$(hostname)}.$$.${RANDOM}
echo "This is the log file for the script." > $LOG_FILE
echo -e "\n pwd is:\n$PWD\n" >> $LOG_FILE
# echo -e "\nScript file is:\n$0\n" >> $LOG_FILE # for regular script usage (no qsub)
echo -e "\nScript file is:\n${JOB_NAME:=$(readlink -m $0)}\n" >> $LOG_FILE 
echo -e "\nScript file contents:\n\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%" >> $LOG_FILE
cat "$0" >> $LOG_FILE
# ~~~~~~~~~~~~ #

# ~~~~~~ # ~~~~~~ # ~~~~~~ run commands ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ 

# ~~~~~~ # ~~~~~~  QUALITY CONTROL # ~~~~~~ # ~~~~~~ # ~~~~~~ 
# # set a QC outdir
tmpQCdir="${OUTDIR}/quality_control"
mkdir -p "$tmpQCdir"

# FastQC; check raw data quality
tmpFastQCdir1="${tmpQCdir}/FastQC_$(basename $FASTQ_R1)"
tmpFastQCdir2="${tmpQCdir}/FastQC_$(basename $FASTQ_R2)"
mkdir -p "$tmpFastQCdir1"
mkdir -p "$tmpFastQCdir2"
fastqc --threads "$THREADS" --nogroup --outdir "$tmpFastQCdir1" "$FASTQ_R1"
fastqc --threads "$THREADS" --nogroup --outdir "$tmpFastQCdir2" "$FASTQ_R2"

# Kraken; check for potential contaminants
tmpKrakendir="${tmpQCdir}/Kraken"
mkdir -p "$tmpKrakendir"
tmp_R1="${tmpKrakendir}/$(basename $FASTQ_R1)"
tmp_R2="${tmpKrakendir}/$(basename $FASTQ_R2)"
# # only use the first 1,000,000 reads !!
zcat $FASTQ_R1 | head -4000000 | gzip > "$tmp_R1"
zcat $FASTQ_R2 | head -4000000 | gzip > "$tmp_R2"
$HOME/software/bin/kraken --gzip-compressed --fastq-input --paired $tmp_R1 $tmp_R2 | $HOME/software/bin/kraken-report | awk -F $'\t' '$1>0.1' > "${tmpKrakendir}/kraken_paired.$(basename $FASTQ_R1).txt"




# ~~~~~~ # ~~~~~~ # ~~~~~~  TRIMMING FASTQ # ~~~~~~ # ~~~~~~ # ~~~~~~ 
# http://hannonlab.cshl.edu/fastx_toolkit/commandline.html
# # trim back the last 1bp on all samples
# need to unzip the fastq for the trimmer, then save a zipped file

# set some local files; raw data is symlinked, some steps are easier to work from pwd instead
tmpFASTQ_R1="${OUTDIR}/$(basename $FASTQ_R1)"
tmpFASTQ_R2="${OUTDIR}/$(basename $FASTQ_R2)"
# local trimmed fastq.gz outputs
trim_R1="${tmpFASTQ_R1//.fastq.gz/}.trim1.fastq.gz"
trim_R2="${tmpFASTQ_R2//.fastq.gz/}.trim1.fastq.gz"

if [ ! -f $trim_R1 ]; then
  echo -e "Trimming reads file \n ${FASTQ_R1}"
  # unzip, trim, re-zip
  gunzip -c "$FASTQ_R1" | $HOME/software/bin/fastx_trimmer -Q33 -t 1 | gzip > "$trim_R1"
fi

if [ ! -f $trim_R2 ]; then
  echo -e "Trimming reads file \n ${FASTQ_R2}"
  gunzip -c "$FASTQ_R2" | $HOME/software/bin/fastx_trimmer -Q33 -t 1 | gzip > "$trim_R2"
fi

# ~~~~~~ # ~~~~~~ # ~~~~~~  ALIGNMENT # ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ 
# align with bowtie2
# /local/data/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome
tmpGenome="/local/data/iGenomes/Homo_sapiens/UCSC/${GENOME}/Sequence/Bowtie2Index/genome"
if [ ! -f ${SAMPLEID}.bam ]; then
  # align with bowtie2, pipe SAM output to samtools for BAM conversion & sorting
  # tee the stderr to save a separate copy of the alignment stats
  (bowtie2 --threads "$THREADS" --local -x "$tmpGenome" -q -1 "$trim_R1" -2 "$trim_R2" | samtools view -@ "$THREADS" -Sb1 - | samtools sort -m 10G -@ "$THREADS" - "$SAMPLEID" ) 3>&1 1>&2 2>&3 | tee ${SAMPLEID}_align-stats.txt 
fi
# for single read use this:
# (bowtie2 --threads "$THREADS" --local -x "$tmpGenome" -q -U "$FASTQ_R1" | samtools view -@ "$THREADS" -Sb1 - | samtools sort -m 10G -@ "$THREADS" - "$SAMPLEID" ) 3>&1 1>&2 2>&3 | tee ${SAMPLEID}_align-stats.txt 

# index the alignment file
if [ ! -f ${SAMPLEID}.bam.bai ]; then
samtools index "${SAMPLEID}.bam"
fi

# get the mapq values for R (make sure this is done before mapq filtering)
if [ ! -f ${SAMPLEID}_mapq-scores_all.txt ]; then
samtools view "${SAMPLEID}.bam" | cut -f5 | sort -n > ${SAMPLEID}_mapq-scores_all.txt
fi

# also make a frequency table for easy text viewing
if [ ! -f ${SAMPLEID}_mapq-scores_freq.txt ]; then
samtools view "${SAMPLEID}.bam" | cut -f5 | sort -n | uniq -c > ${SAMPLEID}_mapq-scores_freq.txt
fi

# make histogram of the mapq frequencies
# in the future update to make sure that the arg passed was a file that ends in .Rnw
Rscript --slave --no-save --no-restore - "${SAMPLEID}_mapq-scores_all.txt" <<EOF
  ## R code
  cat("\nR loaded\n")
  args <- commandArgs(TRUE)
  cat("Script args are:\n")
  args
  pdf("alignment_stats_mapq_hist.pdf",width = 8,height = 8)
  hist(as.numeric(scan(file = args[1],character(0), sep = "\n")),main="MAPQ Alignment Quality Score Distribution",xlab="Score",col="gray")
  dev.off()
EOF
#,labels = TRUE
# btw here are some other ways to pass heredoc to R
# R --slave <<EOF
# Rscript --slave --vanilla - <<EOF






# some user manuals that are useful for these tools:
# http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
# http://www.htslib.org/doc/samtools-0.1.19.html
# http://samtools.github.io/hts-specs/SAMv1.pdf
# https://edwards.sdsu.edu/research/creating-indexed-bam-files-from-bowtie-alignments/




