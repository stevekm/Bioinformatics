#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
input_file <- args[1]
output_file <- args[2]

# read in file
variant_df <- read.delim(file = input_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
save.image("loaded.Rdata")
if (nrow(variant_df) > 0) {
    # split off tumor and normal allelic depth columns
    normal_ad_df <- as.data.frame(do.call('rbind', strsplit(x = variant_df[["NORMAL.AD"]], split = ',')))
    names(normal_ad_df) <- c("NORMAL.AD.ref", "NORMAL.AD.alt")
    normal_ad_df[["NORMAL.AD.total"]] <- as.numeric(as.character(normal_ad_df[["NORMAL.AD.ref"]])) + as.numeric(as.character(normal_ad_df[["NORMAL.AD.alt"]]))
    
    tumor_ad_df <- as.data.frame(do.call('rbind', strsplit(x = variant_df[["TUMOR.AD"]], split = ',')))
    names(tumor_ad_df) <- c("TUMOR.AD.ref", "TUMOR.AD.alt")
    tumor_ad_df[["TUMOR.AD.total"]] <- as.numeric(as.character(tumor_ad_df[["TUMOR.AD.ref"]])) + as.numeric(as.character(tumor_ad_df[["TUMOR.AD.alt"]]))
    
    # add the new columns to the dataset
    new_df <- cbind(variant_df, normal_ad_df)
    new_df <- cbind(new_df, tumor_ad_df)
    
} else {
    new_cols <- c("NORMAL.AD.ref", "NORMAL.AD.alt", "NORMAL.AD.total", "TUMOR.AD.ref", "TUMOR.AD.alt", "TUMOR.AD.total")
    # make an empty dataframe with new columns
    d2 <- do.call('cbind', lapply(X = new_cols, FUN = function(x){
        d <- as.data.frame(character())
        names(d) <- x
        return(d)
    }))
    new_df <- cbind(variant_df, d2)
}


# save output file
write.table(x = new_df, file = output_file, sep = '\t', row.names = FALSE, col.names = TRUE)

save.image("finish.Rdata")