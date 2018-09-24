#!/usr/bin/env Rscript
args = commandArgs(TRUE)
input_maf = args[1]
output_maf = args[2]

filter_TCGAmaf <- function(df){
    # inspired by:
    # https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-017-0424-2
    # https://doi.org/10.1186/s13073-017-0424-2
    # Analysis of 100,000 human cancer genomes reveals the landscape of tumor mutational burden. Chalmers ZR et. al.
    
    # keep entries not listed as non-coding, and have a protein coding entry
    df <- droplevels(df[ grep(pattern = '*non_coding*', x = df[["Consequence"]], invert = TRUE), ])
    df <- droplevels(df[ which(df[["HGVSp"]] != ''), ])
    
    # keep entries that lack COSMIC identifers
    df <- droplevels(df[ which(df[["COSMIC"]] == ''), ])
    
    # keep entries that lack dbSNP identifiers or are listed as 'novel'
    df <- droplevels(df[ which(df[["dbSNP_RS"]] == '' | df[["dbSNP_RS"]] == 'novel'), ])
    
    # only keep entries with 'NA' in ExAC fields
    df <- droplevels(df[ is.na(df[["ExAC_AF"]]), ])

    # add a column for the coding change notation for later
    df[["AAChange"]] <- sprintf('%s>%s', df[["Reference_Allele"]], df[["Tumor_Seq_Allele1"]])
    
    return(df)
}

df <- read.delim(file = input_maf, header = TRUE, sep = '\t', comment.char = '#')

df <- filter_TCGAmaf(df)

write.table(x = df, file = output_maf, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)

save.image()