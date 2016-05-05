#!/bin/bash
# print the script commands to stderr on execution; stderr gets saved in the qsub log files
set -x
##
## USAGE: PARCLIP-wavClusteR-pipeline_singleRead.sh /path/to/outdir /path/to/input_file.fastq.gz <sampleID> <ref_genome> /path/to/ref_genome  
## This script will run PAR-CLIP analysis pipeline with bowtie paired end alignment and wavClusteR analysis pipeline
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
# these need to be installed
module load bowtie/1.0.0
module load samtools/1.2.1
module unload r
module load r/3.2.3 

# get the script arguments
OUTDIR="$1" # outdir
FASTQ_R1="$2" # Read 1 input file
# FASTQ_R2="$3" # Read 2 file
SAMPLEID="$3" # ID for the sample... 
GENOME="$4" # the name of the genome to be used e.g. hg19, mm10, etc
GENOME_PATH="$5" # the path to the reference genome to use for Bowtie; e.g. /local/data/iGenomes/Mus_musculus/UCSC/mm10/Sequence/BowtieIndex/genome
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

# exit

# ~~~~~~ # ~~~~~~ # ~~~~~~ run commands ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ 
# our files are fastq.gz, so need to unzip them for this... :(
# if [ ! -f ${SAMPLEID}.1.fastq ]; then zcat "$FASTQ_R1" > "${SAMPLEID}.1.fastq"; fi
# if [ ! -f ${SAMPLEID}.2.fastq ]; then zcat "$FASTQ_R2" > "${SAMPLEID}.2.fastq"; fi

# align with bowtie
if [ ! -f ${SAMPLEID}.sam ]; then 
# (bowtie "$GENOME_PATH" --threads "$THREADS" -v 2 -m 10 --best --strata -1 "${SAMPLEID}.1.fastq" -2 "${SAMPLEID}.2.fastq" -S "${SAMPLEID}.sam" ) 2>${SAMPLEID}.sam.align_stats.txt
# bowtie "$GENOME_PATH" --threads "$THREADS" -v 2 -m 10 --best --strata -1 "${SAMPLEID}.1.fastq" -2 "${SAMPLEID}.2.fastq" -S "${SAMPLEID}.sam" 2>${SAMPLEID}.sam.align_stats.txt
(zcat ${FASTQ_R1} | bowtie "$GENOME_PATH" --threads "$THREADS" -v 2 -m 10 --best --strata - -S "${SAMPLEID}.sam" ) 2>${SAMPLEID}.sam.align_stats.txt
fi

# convert SAM to BAM, sort
if [ ! -f ${SAMPLEID}.bam ]; then
cat ${SAMPLEID}.sam | samtools view -@ "$THREADS" -Sb1 - | samtools sort -m 10G -@ "4" - "${SAMPLEID}"
fi

# index the BAM
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
if [ ! -f alignment_stats_mapq_hist ]; then
Rscript --slave --no-save --no-restore - "${SAMPLEID}_mapq-scores_all.txt" <<EOFA
  ## R code
  cat("\nR loaded\n")
  args <- commandArgs(TRUE)
  cat("Script args are:\n")
  args
  pdf("alignment_stats_mapq_hist.pdf",width = 8,height = 8)
  hist(as.numeric(scan(file = args[1],character(0), sep = "\n")),main="MAPQ Alignment Quality Score Distribution",xlab="Score",col="gray")
  dev.off()
EOFA
#,labels = TRUE # btw here are some other ways to pass heredoc to R # R --slave <<EOF # Rscript --slave --vanilla - <<EOF
fi

# run the wavClusteR pipeline, in R
(
Rscript --slave --no-save --no-restore - "${SAMPLEID}" "${SAMPLEID}.bam" <<EOFB
  ## R code
  cat("\nR loaded\n")
  args <- commandArgs(TRUE)
  cat("Script args are:\n")
  args
  
  sampleID<-args[1]
  filename<-args[2]
  
  # load wavClusteR package
  cat("\nLoading wavClusteR package\n")
  library(wavClusteR)
  
  # load the bam file
  cat("\nReading in the sorted bam file\n")
  Bam <- readSortedBam(filename = filename)
  
  # Extracting informative genomic positions
  cat("\nExtracting informative genomic positions\n")
  countTable <- getAllSub( Bam, minCov = 10 )
  
  #  Inspecting the substitutions profile (diagnostic I)
  cat("\nInspecting the substitutions profile\n")
  cat("\nSaving plot ",paste0(sampleID,".count_table.pdf"),"\n")
  pdf(paste0(sampleID,".count_table.pdf"),width = 9,height = 9)
  plotSubstitutions( countTable, highlight = "TC" )
  dev.off()
  
  # Estimating the model
  cat("\nEstimating the model\n")
  model <- fitMixtureModel(countTable, substitution = "TC")
  data(model)
  str(model)
  
  cat("\nSaving plot ",paste0(sampleID,".model_estimate.pdf"),"\n")
  pdf(paste0(sampleID,".model_estimate.pdf"),width = 9,height = 9)
  support <- getExpInterval( model, bayes = TRUE )
  dev.off()
  
  cat("\nSaving plot ",paste0(sampleID,".model_estimate_prob0.9.pdf"),"\n")
  pdf(paste0(sampleID,".model_estimate_prob0.9.pdf"),width = 9,height = 9)
  support <- getExpInterval( model, bayes = FALSE, leftProb = 0.9, rightProb =   0.9 )
  dev.off()
  
  cat("\nSaving plot ",paste0(sampleID,".count_table.pdf"),"\n")
  pdf(paste0(sampleID,".count_table.pdf"),width = 9,height = 9)
  plotSubstitutions( countTable, highlight = "TC", model )
  dev.off()
  
  # Selecting high-confidence PAR-CLIP induced transitions
  cat("\nSelecting high-confidence PAR-CLIP induced transitions\n")
  highConfSub <- getHighConfSub( countTable, 
                               support = support, 
                               substitution = "TC" )
  head( highConfSub )
  
  coverage <- coverage( Bam )
  coverage$chrX
  
  clusters <- getClusters( highConfSub = highConfSub,
                         coverage = coverage,
                         sortedBam = Bam,
                         method = "mrn",
                         threshold = 1,
                         cores = 6 )
  clusters
  
  
  # Merging clusters
  cat("\nMerging clusters\n")
  
  # source("http://bioconductor.org/biocLite.R")
  # biocLite("BSgenome.Hsapiens.UCSC.hg19")
  # ~/R/x86_64-pc-linux-gnu-library/3.2
  # library(BSgenome.Hsapiens.UCSC.hg19)
  # source("http://bioconductor.org/biocLite.R")
  # biocLite("BSgenome.Mmusculus.UCSC.mm10")
  
  cat("\nLoading UCSC\n")
  library(BSgenome.Mmusculus.UCSC.mm10)
  # genome <- BSgenome.Mmusculus.UCSC.mm10
  
  wavclusters <- filterClusters( clusters = clusters, 
                               highConfSub = highConfSub,
                               coverage = coverage, 
                               model = model, 
                               genome = Mmusculus, 
                               # genome = genome, 
                               # genome = Hsapiens, # Mmusculus
                               refBase = "T", 
                               minWidth = 12)
  wavclusters
  
  
  # Post-processing of the binding sites
  # 3.9.1 Exporting substitutions, wavClusters and coverage function
  cat("\nPost-processing of the binding sites\n")
  cat("\nExporting substitutions, wavClusters and coverage function\n")
  exportHighConfSub( highConfSub = highConfSub, 
                   filename = paste0(sampleID,".hcTC.bed"), 
                   trackname = "hcTC", 
                   description = "hcTC" )
  exportClusters( clusters = wavclusters, 
                filename = paste0(sampleID,".wavClusters.bed"), 
                trackname = "wavClusters", 
                description = "wavClusters" )
  exportCoverage( coverage = coverage, filename = paste0(sampleID,".coverage.bigWig" ))
  
  
  #  Annotating binding sites
  # source("https://bioconductor.org/biocLite.R")
  # biocLite("GenomicFeatures")
  cat("\nAnnotating binding sites\n")
  cat("\nLoading GenomicFeatures package\n")
  library("GenomicFeatures")
  # txDB <- makeTxDbFromUCSC(genome = "hg19", tablename = "ensGene")
  txDB <- makeTxDbFromUCSC(genome = "mm10", tablename = "ensGene")
  
  #cat("\nSaving plot ",paste0(sampleID,".annotate_clusters.pdf"),"\n")
  #pdf(paste0(sampleID,".annotate_clusters.pdf"),width = 9,height = 9)
  #annotateClusters( clusters = wavclusters, 
  #                txDB = txDB, 
  #                plot = TRUE, 
  #                verbose = TRUE)
  #dev.off()
  
  # Computing cluster metagene profiles
  cat("\nComputing cluster metagene profiles\n")
  cat("\nSaving plot ",paste0(sampleID,".getMetaGene.pdf"),"\n")
  pdf(paste0(sampleID,".getMetaGene.pdf"),width = 9,height = 9)
  getMetaGene( clusters = wavclusters, 
             txDB = txDB, 
             upstream = 1e3, 
             downstream = 1e3, 
             nBins = 40, 
             nBinsUD = 10, 
             minLength = 1, 
             plot = TRUE,
             verbose = TRUE ) 
  dev.off()
  
  cat("\nSaving plot ",paste0(sampleID,".getMetaTSS.pdf"),"\n")
  pdf(paste0(sampleID,".getMetaTSS.pdf"),width = 9,height = 9)
  getMetaTSS( sortedBam = Bam, 
            txDB = txDB, 
            upstream = 1e3, 
            downstream = 1e3, 
            nBins = 40, 
            unique = FALSE, 
            plot = TRUE, 
            verbose = TRUE ) 
  dev.off()
  
  # Computing the size distribution of wavClusters
  cat("\nComputing the size distribution of wavClusters\n")
  cat("\nSaving plot ",paste0(sampleID,".size_dist_scatter.pdf"),"\n")
  pdf(paste0(sampleID,".size_dist_scatter.pdf"),width = 9,height = 9)
  plotSizeDistribution( clusters = wavclusters, showCov = TRUE, col = "skyblue2" )
  dev.off()
  
  cat("\nSaving plot ",paste0(sampleID,".size_dist_hist.pdf"),"\n")
  pdf(paste0(sampleID,".size_dist_hist.pdf"),width = 9,height = 9)
  plotSizeDistribution( clusters = wavclusters, col = "skyblue2" )
  dev.off()
  
  # Visualizing wavCluster statistics and meta data
  cat("\nVisualizing wavCluster statistics and meta data\n")
  cat("\nSaving plot ",paste0(sampleID,".cluster_stats.pdf"),"\n")
  pdf(paste0(sampleID,".cluster_stats.pdf"),width = 9,height = 9)
  plotStatistics( clusters = wavclusters, 
                corMethod = "spearman", 
                lower = panel.smooth )
  dev.off()
  
  # sink("R_sessionInfo.txt")
  system('uname -srv',intern=T)
  sessionInfo()
  # sink()
  save.image(compress = TRUE)
EOFB
) 1>&2 # write the stdout to stderr




