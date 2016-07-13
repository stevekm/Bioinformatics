#!/bin/bash
# print the script commands to stderr on execution; stderr gets saved in the qsub log files
set -x
##
## USAGE: analysis_pipeline.sh /path/to/outdir /path/to/input_file_R1.fastq.gz /path/to/input_file_R2.fastq.gz <sampleID> <ref_genome> /path/to/ref_genome  
## This script will run the analysis pipeline; submit this script with qsub!
## 

# ~~~~~~ permissions ~~~~~~ #
# To make stuff group-writeable
umask 007

# ~~~~~~ script args processing ~~~~~~ #
# # if not enough args, output USAGE info lines (which start with '##') and exit
if (($# != 6)); then
  grep '^##' $0
  exit
fi

# ~~~~~~ params ~~~~~~ #
# put programs to load here, example:
# module load bowtie/1.0.0
# module load samtools/1.2.1
module unload r
module load r/3.2.3 

# ~~~~~~ script args ~~~~~~ #
# get the script arguments here
OUTDIR="$1" # outdir
FASTQ_R1="$2" # Read 1 input file
FASTQ_R2="$3" # Read 2 file
SAMPLEID="$4" # ID for the sample... 
GENOME="$5" # the name of the genome to be used e.g. hg19, mm10, etc
GENOME_PATH="$6" # the path to the reference genome to use for Bowtie; e.g. /local/data/iGenomes/Mus_musculus/UCSC/mm10/Sequence/BowtieIndex/genome
# THREADS is defined by qsub; if not, set with 8 threads
THREADS=${NSLOTS:=8} 
echo "OUTDIR is $OUTDIR"
echo "FASTQ_R1 is $FASTQ_R1"
echo "FASTQ_R2 is $FASTQ_R2"
echo "SAMPLEID is $SAMPLEID"
echo "GENOME is $GENOME"
echo "GENOME_PATH is $GENOME_PATH"
echo "THREADS is $THREADS"

# ~~~~~~ Script Logging ~~~~~~~ #
# save a hardcopy of this script file 
zz=$(basename $0)
# get the script dir
za=$(dirname $0)
# set the script log file, use JOBNAME if its set by qsub
if [ ! -d logs ]; then mkdir logs; fi
LOG_FILE=logs/scriptlog.${JOB_NAME:=$(basename $0)}.$(date -u +%Y%m%dt%H%M%S).${HOSTNAME:-$(hostname)}.$$.${RANDOM}
echo "This is the log file for the script." > $LOG_FILE
echo -e "\n pwd is:\n$PWD\n" >> $LOG_FILE
# echo -e "\nScript file is:\n$0\n" >> $LOG_FILE # for regular script usage (no qsub)
echo -e "\nScript file is:\n${JOB_NAME:=$(readlink -m $0)}\n" >> $LOG_FILE 
echo -e "\nScript file contents:\n\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%" >> $LOG_FILE
cat "$0" >> $LOG_FILE
# ~~~~~~~~~~~~ #


# ~~~~~~ # ~~~~~~ # ~~~~~~ run commands ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ 

# place your analysis commands to run on each sample here

# make some sample plots, from here: http://www.r-bloggers.com/normal-distribution-functions/
(
Rscript --slave --no-save --no-restore - "${SAMPLEID}" <<EOFB
  ## R code
  cat("\nR loaded\n")
  # get R script args
  args <- commandArgs(TRUE)
  cat("Script args are:\n")
  args
  
  # get Sample ID
  sampleID<-args[1]
  
  # make some sample plots
  # # set some random numbers
  set.seed(3000)
  xseq<-seq(-4,4,.01)
  densities<-dnorm(xseq, 0,1)
  cumulative<-pnorm(xseq, 0, 1)
  randomdeviates<-rnorm(1000,0,1)
  
  # write out some stats
  sink("R_stats.txt")
  cat(paste0(sampleID),"\n")
  cat("Some sample stats\n")
  head(xseq)
  head(densities)
  head(cumulative)
  head(randomdeviates)
  sink()
  
  # save the first plot
  cat("\nSaving plot ",paste0(sampleID,".random_distribution.pdf"),"\n")
  pdf(paste0(sampleID,".random_distribution.pdf"),width = 9,height = 9)
  par(mfrow=c(1,3), mar=c(3,4,4,2))
  plot(xseq, densities, col="darkgreen",xlab="", ylab="Density", type="l",lwd=2, cex=2, main="PDF of Standard Normal", cex.axis=.8)
  plot(xseq, cumulative, col="darkorange", xlab="", ylab="Cumulative Probability",type="l",lwd=2, cex=2, main="CDF of Standard Normal", cex.axis=.8)
  hist(randomdeviates, main="Random draws from Std Normal", cex.axis=.8, xlim=c(-4,4))
  dev.off()
  
  # make some more random numbers
  xseq<-seq(-4,4,.01)
  y<-2*xseq + rnorm(length(xseq),0,5.5)
  
  # make another plot
  cat("\nSaving plot ",paste0(sampleID,".distribution_histogram.pdf"),"\n")
  pdf(paste0(sampleID,".distribution_histogram.pdf"),width = 9,height = 9)
  hist(y, prob=TRUE, ylim=c(0,.06), breaks=20)
  curve(dnorm(x, mean(y), sd(y)), add=TRUE, col="darkblue", lwd=2)
  dev.off()
  
  # sink("R_sessionInfo.txt")
  system('uname -srv',intern=T)
  sessionInfo()
  # sink()
  save.image(compress = TRUE)
  
EOFB
) 1>&2 # write the stdout to stderr

