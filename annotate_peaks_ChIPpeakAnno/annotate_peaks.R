#!/usr/bin/env Rscript

## USAGE: annotate_peaks.R /path/to/peaks.bed /path/to/output/annotated_peaks.tsv
## DESCRIPTION: This script will run annotate peaks with ChIPpeakAnno, using hg19

# get script args
args <- commandArgs(TRUE)

cat("\nScript args are:\n")
print(args)

input_peaks_file <- args[1]
output_annotated_peaks_file <- args[2]

cat("\nLoading packages...\n")

# source("https://bioconductor.org/biocLite.R")
# biocLite("ChIPpeakAnno")
library(ChIPpeakAnno)
library(biomaRt)

# read in the BED file
cat("\nReading in the BED file...\n")
peaks_granges <- toGRanges(input_peaks_file, format="BED", header=FALSE) 

# for hg19
# get biomart reference genome information
# cat("\nGetting biomart reference information for hg19...\n")
# martEns <- useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", verbose=F)
# martEnsTSS <- getAnnotation(mart=martEns, featureType="TSS")
# martEnsDF <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "gene_biotype"), mart=martEns)
# save(martEns, martEnsTSS, martEnsDF, file = "biomart_data.RData")

cat("\nLoading biomart reference information for hg19...\n")
load("biomart_data.RData")

# get the annotations
cat("\nGetting annotations...\n")
peaks_granges <- annotatePeakInBatch(peaks_granges, AnnotationData = martEnsTSS, PeakLocForDistance = "middle", FeatureLocForDistance = "TSS", output = "shortestDistance", multiple = TRUE)

# merge the annotations with the peaks
cat("\nMerging annotations...\n")
peaks_granges_df <- merge(as.data.frame(peaks_granges) , martEnsDF , by.x=c("feature"), by.y=c("ensembl_gene_id") , all.x=TRUE)

# !! NOTE: Need to subtract 1 from the 'start' to get the original coordinate !!

# save the output
cat("\nSaving the output...\n")
write.table(peaks_granges_df, row.names = FALSE, sep = '\t', quote = FALSE, 
            file = output_annotated_peaks_file)


sessionInfo()
