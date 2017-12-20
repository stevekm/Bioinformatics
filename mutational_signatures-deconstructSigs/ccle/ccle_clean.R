#!/usr/bin/env Rscript

# script to clean the Broad CCLE public RNASeq tumor dataset from here:
# https://portals.broadinstitute.org/ccle/data
# for use with "deconstructSigs" to make downstream scripts easier

# Run this first!

library("BSgenome.Hsapiens.UCSC.hg19")
library("deconstructSigs")


# ~~~~~ LOAD DATA ~~~~~ # 
ccle_data_file <- "ccle2maf_081117.txt"
ccle <- read.delim(file = ccle_data_file, header = TRUE, sep = '\t', check.names = FALSE)


# ~~~~~ CONVERT DATA ~~~~~ # 
# fix the chrom names
ccle[["Chromosome"]] <- sprintf('chr%s', as.character(ccle[["Chromosome"]]))

# keep only entries with chroms in the reference data
ccle <- ccle[which(ccle[["Chromosome"]] %in% seqnames(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)), ]

# need at least 50 variants per sample, remove samples that dont have enough variants
# actually looks like you need at least 55 variants? 
tumor_samples <- data.frame(table(ccle[["Tumor_Sample_Barcode"]]))
remove_samples <- as.character(tumor_samples[which(tumor_samples[["Freq"]] < 55), "Var1"])
ccle <- ccle[which(! as.character(ccle[["Tumor_Sample_Barcode"]]) %in% remove_samples), ]
ccle <- droplevels(ccle)

save(ccle, file = "ccle_clean.Rdata", compress = TRUE)


# convert to signatures format
sigs.input <- mut.to.sigs.input(mut.ref = ccle, 
                                sample.id = "Tumor_Sample_Barcode", 
                                chr = "Chromosome", 
                                pos = "Start_position", 
                                ref = "Reference_Allele", 
                                alt = "Tumor_Seq_Allele1")


# save the final object & table
save(sigs.input, file = "ccle_sigs.input.Rdata", compress = TRUE)
write.table(x = sigs.input, file = "ccle_sigs.input.tsv", quote = FALSE, sep = '\t', row.names = TRUE, col.names = TRUE)

# save just the sample IDs
fileConn <- file("ccle_samples.txt")    
writeLines(rownames(x = sigs.input), fileConn)    
close(fileConn)

# read the samples from the text file
# fileConn <- file("ccle_samples.txt")    
# sampleIDs <- readLines(fileConn)
# close(fileConn)
