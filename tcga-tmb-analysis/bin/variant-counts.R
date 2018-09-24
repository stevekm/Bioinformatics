#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
input_maf <- args[1]
output_counts <- args[2]
output_aachange <- args[3]
study <- args[4]
project <- args[5]
caller <- args[6]

# load df
df <- read.delim(input_maf, sep = "\t")

# add a column to count on
df[["count"]] <- 1

# take the aggregate sum of samples
agg <- aggregate( count ~ Tumor_Sample_Barcode + tumor_bam_uuid, data = df, FUN = sum)

# add some extra columns
agg[["study"]] <- study
agg[["project"]] <- project
agg[["caller"]] <- caller

# save file
write.table(x = agg, file = output_counts, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

# get an aggregate of amino acid changes
aa_agg <- aggregate( count ~ AAChange, data = df, FUN = sum)
aa_agg[["study"]] <- study
aa_agg[["project"]] <- project
aa_agg[["caller"]] <- caller

# only keep the known transition and inversions
aa_to_keep <- c('C>T', 'C>A', 'C>G', 'T>C', 'T>A', 'T>G')
aa_agg <- aa_agg[aa_agg[["AAChange"]] %in% aa_to_keep, ]

write.table(x = aa_agg, file = output_aachange, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

save.image()
