#!/usr/bin/env Rscript

# compare the signatures generated for the CCLE dataset by the previous scripts

# ~~~~~ LOAD DATA ~~~~~ # 
# dir containg the signature .Rds files
signatures_dir <- "signatures"

# named vector of paths to signature files
samples <- setNames(nm = sapply(X = strsplit(x = dir(signatures_dir, 
                                          pattern = "*.Rds"), 
                                  split = '.', 
                                  fixed = TRUE), 
                     FUN = function(x) x[1]), 
         object = dir(signatures_dir, 
                      full.names = TRUE,
                      pattern = "*.Rds"))

# list to hold all the data
sample_signatures <- list()

for(i in seq_along(samples)){
    sampleID <- names(samples)[i]
    signatures_file <- samples[[i]]
    
    # load the .Rds
    message(sprintf("loading: %s\n", signatures_file))
    signatures <- readRDS(file = signatures_file)
    sample_signatures[[sampleID]] <- signatures
}



# ~~~~~ EXTRACT DATA ~~~~~ # 
# just the all the Signature weights from the samples
all_weights <- lapply(sample_signatures, function(x) x[["weights"]])
weights_df <- do.call(rbind, all_weights)

# sort by the desired signature value
weights_df <- weights_df[with(weights_df, order(-Signature.3)), ]

# save to file
write.table(x = weights_df, file = 'ccle_signature_weights.tsv', quote = FALSE, sep = '\t', row.names = TRUE, col.names = NA)
