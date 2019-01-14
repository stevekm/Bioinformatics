#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
input_file <- args[1]
output_file <- args[2]

# read in file
variant_df <- read.delim(file = input_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)

# split off tumor and normal allelic depth columns
normal_ad_df <- as.data.frame(do.call('rbind', strsplit(x = variant_df[["NORMAL.AD"]], split = ',')))
names(normal_ad_df) <- c("NORMAL.AD.ref", "NORMAL.AD.alt")
normal_ad_df[["NORMAL.AD.total"]] <- as.numeric(as.character(normal_ad_df[["NORMAL.AD.ref"]])) + as.numeric(as.character(normal_ad_df[["NORMAL.AD.alt"]]))
normal_ad_df[["NORMAL.AD.total.x5"]] <- normal_ad_df[["NORMAL.AD.total"]] * 5
normal_ad_df[["NORMAL.AD.alt_div_NORMAL.AD.total"]] <- as.numeric(as.character(normal_ad_df[["NORMAL.AD.alt"]])) / normal_ad_df[["NORMAL.AD.total"]]

tumor_ad_df <- as.data.frame(do.call('rbind', strsplit(x = variant_df[["TUMOR.AD"]], split = ',')))
names(tumor_ad_df) <- c("TUMOR.AD.ref", "TUMOR.AD.alt")
tumor_ad_df[["TUMOR.AD.total"]] <- as.numeric(as.character(tumor_ad_df[["TUMOR.AD.ref"]])) + as.numeric(as.character(tumor_ad_df[["TUMOR.AD.alt"]]))
tumor_ad_df[["TUMOR.AD.alt_div_TUMOR.AD.total"]] <- as.numeric(as.character(tumor_ad_df[["TUMOR.AD.alt"]])) / tumor_ad_df[["TUMOR.AD.total"]]

tumor_ad_df[["TUMOR.AD.alt_div_TUMOR.AD.total_gt_NORMAL.AD.alt_div_NORMAL.AD.total"]] <- tumor_ad_df[["TUMOR.AD.alt_div_TUMOR.AD.total"]] > (normal_ad_df[["NORMAL.AD.alt_div_NORMAL.AD.total"]] * 5)

# ( Normal Allelic depth alt / ( Normal Allelic depth ref + Normal Allelic depth alt ) )  < 0.05
# ( Tumor Allelic depth alt / (Tumor Allelic depth ref + Tumor Allelic depth alt ) )  > 0.03
# Tumor Allelic depth alt > 5
# ( Tumor Allelic depth alt / ( Tumor Allelic depth ref + Tumor Allelic depth alt ) ) > ( Normal Allelic depth alt / ( Normal Allelic depth ref + Normal Allelic depth alt ) ) * 5

# add the new columns to the dataset
new_df <- cbind(variant_df, normal_ad_df)
new_df <- cbind(new_df, tumor_ad_df)

# save output file
write.table(x = new_df, file = output_file, sep = '\t', row.names = FALSE, col.names = TRUE)

save.image("loaded.Rdata")