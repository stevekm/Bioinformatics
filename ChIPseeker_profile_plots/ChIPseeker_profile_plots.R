#!/usr/bin/env Rscript

## USAGE: ChIPseeker_profile_plots.R sample_ID /path/to/peaks.bed /path/to/output.pdf 
## DESCRIPTION: This script will create profile plots of ChIP-Seq peaks
## using ChIPSeeker
## thanks to Olivier
## http://docplayer.net/51086377-Chip-sequencing-analysis-in-r-bioconductor.html

# ~~~~~ LOAD PACKAGES ~~~~~ #
library("optparse")
library("ggplot2")
library("ChIPseeker")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")


# ~~~~~ CUSTOM FUNCTIONS ~~~~~ #
make_profile_plot <- function(sampleID, peaks_file, output_file){
    # make the genomic regions profile plot for a peaks file
    
    # read in the peaks
    peaks <- readPeakFile(peaks_file)
    
    # make the matrix
    tagMatrix <- getTagMatrix(peaks, weightCol=NULL, windows=promoter)
    
    # make the plot
    profile_plot <- plotAvgProf(tagMatrix, 
                                xlim=c(-3000, 3000), 
                                xlab="Genomic Region (5'->3')", 
                                ylab = "Read Count Frequency")
    profile_plot <- profile_plot + ggtitle(sampleID)
    
    # save the plot
    pdf(file = output_file, height = 10, width = 10 )
    print(profile_plot)
    dev.off()
}




# ~~~~~ GET SCRIPT ARGS ~~~~~ #
args <- commandArgs(TRUE)

sampleID <- args[1]
peaks_file <- args[2]
output_file <- args[3]



# ~~~~~ RUN ~~~~~ #
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoter <- getPromoters(TxDb = txdb, upstream=3000, downstream=3000)
make_profile_plot(sampleID, peaks_file, output_file)
