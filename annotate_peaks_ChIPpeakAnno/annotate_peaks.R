#!/usr/bin/env Rscript

## USAGE: annotate_peaks.R /path/to/peaks.bed /path/to/output/annotated_peaks.tsv
## DESCRIPTION: This script will run annotate peaks with ChIPpeakAnno, using hg19

# get script args
args <- commandArgs(TRUE)

message("\nScript args are:\n")
print(args)

input_peaks_file <- args[1]
output_annotated_peaks_file <- args[2]

message("\nLoading packages...\n")

# source("https://bioconductor.org/biocLite.R")
# biocLite("ChIPpeakAnno")
library(ChIPpeakAnno)
library(biomaRt)

# read in the BED file
message("\nReading in the BED file...\n")
peaks_granges <- toGRanges(input_peaks_file, format="BED", header=FALSE) 

# for hg19
# get biomart reference genome information
# check for a saved copy first..
biomart_data_file <- file.path(getwd(), "biomart_data.RData")
if(file.exists(biomart_data_file)){
    message(sprintf("Found biomaRt data file:\n%s\nLoading data from file...", biomart_data_file))
    load(biomart_data_file)
} else {
    message("Saved biomaRt data file not found!")
    message("Retreiving reference information for hg19 from biomaRt, this might take a few minutes...")
    martEns <- useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", verbose=F)
    martEnsTSS <- getAnnotation(mart=martEns, featureType="TSS")
    martEnsDF <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "gene_biotype"), mart=martEns)
    message(sprintf("Saving biomaRt data to file:\n%s\n", biomart_data_file))
    save(martEns, martEnsTSS, martEnsDF, file = biomart_data_file)
}


# get the annotations
message("\nGetting annotations...\n")
peaks_granges <- annotatePeakInBatch(peaks_granges, AnnotationData = martEnsTSS, PeakLocForDistance = "middle", FeatureLocForDistance = "TSS", output = "shortestDistance", multiple = TRUE)

# merge the annotations with the peaks
message("\nMerging annotations...\n")
peaks_granges_df <- merge(as.data.frame(peaks_granges) , martEnsDF , by.x=c("feature"), by.y=c("ensembl_gene_id") , all.x=TRUE)

# save the output
message("\nSaving the output...\n")
write.table(peaks_granges_df, row.names = FALSE, sep = '\t', quote = FALSE, 
            file = output_annotated_peaks_file)


message("Session Information:\n")
sessionInfo()
