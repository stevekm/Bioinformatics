#!/usr/bin/env Rscript
args <- commandArgs(T)

signaturesRds <- args[1]
output_tsv <- args[2]

signatures <- readRDS(signaturesRds)
save.image("loaded.Rdata")

weights <- signatures[["weights"]]

write.table(x = weights, file = output_tsv, sep = '\t', col.names = TRUE, row.names = FALSE)

save.image("finished.Rdata")