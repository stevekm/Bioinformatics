#!/usr/bin/env Rscript

# !! second attempt at this pipeline to try and make things clearer from the start !!

## USAGE: diffbind_geneExpression_analysis.R /path/to/outdir /path/to/microarray_gene_expression.tsv /path/to/diffbind.csv <HistMark> <analysis_params_branch>

## OUTPUT: <outdir>/<SampleID>/<HistMark>/<analysis_params_branch>/[paired_boxplot.pdf genes_up.tsv genes_down.tsv]
## e.g. <outdir>/ZNK/INPUT/peaks.by_sample.macs_broad-cutoff0025

## DESCRIPTION:
# This pipeline script will compare microarray gene expression data with
# DiffBind differential peak binding data
# for each sample, the script will find the set of genes that are up or down regulated in the microarray data
# and create subsets of genes with differentially bound peaks from the DiffBind data
# and plot the DiffBind values for both subsets (up and down regulated genes) in a paired beanplot
# the subset diffbind values will be output as TSV, and the plot


# set up the analysis_output lke this:
# ZNK/INPUT/peaks.by_sample.macs_broad-cutoff0025
# ~~~~~~~~~~~~~~~ #


# ~~~~~~ GET SCRIPT ARGS ~~~~~~~~ #
args <- commandArgs(TRUE); cat("Script args are:\n"); args

# location for the output
main_outdir <- args[1]
setwd(main_outdir)

# path to the TSV file containing the microarray gene expression data
gene_expr_file <- args[2]

# path to the file containing the differential binding DiffBind sheet
diffbind_file <- args[3]

# get the unique ID to use for the sample
histone_mark <- args[4]

# get the analysis params branch for the samples
params_branch <- args[5]
# ~~~~~~~~~~~~~~~ #



# ~~~~~~ READ IN FILE ~~~~~~~~ #
# diffbind data
diffbind_df <- read.csv(diffbind_file)
# gene expression data
gene_expr_df <- read.table(gene_expr_file,sep = '\t',header = TRUE,row.names = 1,quote = "")
# ~~~~~~~~~~~~~~~ #


# ~~~~~ GET PATIENT IDS ~~~~~~~~ #
# FROM GENE EXPRESSION DATA
# get colnames from gene expr table; need to parse out patient names
gene_expr_patientIDs <- colnames(gene_expr_df)
# regex = starts with any number of alpha followed by any number of punct then ends with any number of anything; keep only the first alpha match, unique entries
# # "AGK.Exp.D" - > "AGK"
gene_expr_patientIDs <- unique(gsub(perl = TRUE,pattern = "^([[:alpha:]]*)[[:punct:]]*.*$",replacement = "\\1",x = gene_expr_patientIDs))
#  [1] "AGK" "BVI" "CBK" ...

# MATCH TO DIFFBIND DATA
# figure out which colnames are also present in the diffbind table
# # concatenate the gene expr ID's into search pattern regex to use in the diffbind 
search_pattern <- paste(gene_expr_patientIDs,collapse = "|")
# [1] "AGK|BVI|CBK ... 

# # search for matching colnames
colname_matches <- grep(pattern = search_pattern,x = colnames(diffbind_df),value = TRUE)
# [1] "AGK.D.H3K27AC" "BVI.D.H3K27AC" "CBK.D.H3K27AC" ... 

# # get just the patient ID's back out from DiffBind
patientID_matches <- unique(gsub(perl = TRUE,pattern = "^([[:alpha:]]*)[[:punct:]]*.*$",replacement = "\\1",x = colname_matches))
# [1] "AGK" "BVI" "CBK"
# ~~~~~~~~~~~~~~~ #


# ~~~~~ GET GENE NAMES ~~~~~~~~ #
# get gene expression gene names
gene_expr_gene_names <- rownames(gene_expr_df)
# get genes from Diffbind
diffbind_gene_names <- unique(diffbind_df[["external_gene_name"]])
# grep(x = gene_expr_gene_names,pattern = "^NA",value = TRUE) # check for NA's.. 
# ~~~~~~~~~~~~~~~ #



# ~~~~~ GET COMMON GENES ~~~~~~~~ #
# get the genes in common between the diffbind and gene expression data
gene_names_overlap <- intersect(gene_expr_gene_names,diffbind_gene_names)

# subset the diffbind and gene_expr_df for the entries that match gene expression gene names
diffbind_df_intersect <- droplevels(diffbind_df[diffbind_df[["external_gene_name"]] %in% gene_names_overlap,])
gene_expr_df_intersect <- droplevels(gene_expr_df[rownames(gene_expr_df)  %in% gene_names_overlap,])
# ~~~~~~~~~~~~~~~ #


# ~~~~~~ GET THE CLOSEST DIFFBIND PEAK PER GENE ~~~~~~~~~ #
# in Diffbind subset, getthe peaks closest to each gene
# # for each gene name, find the peak with the lowest "shortestDistance" value
library(data.table)
# make it a data.table type
setDT(diffbind_df_intersect)
# get min value per factor, also remove duplicates
diffbind_df_intersect_min <- diffbind_df_intersect[, .SD[which.min(shortestDistance)], by=external_gene_name]
# diffbind_df_intersect_min <- diffbind_df_intersect[, .SD[shortestDistance == min(shortestDistance)], .(external_gene_name)]
# make sure there are no NA's..
# grep(x=diffbind_df_intersect_min[["external_gene_name"]],pattern = "^NA",value = TRUE)

# remove extraneous columns
diffbind_cols_to_remove <-c("seqnames","start","end","width","strand","Conc","Conc_D","Conc_R","Fold",
                            "p.value","FDR","feature","gene_biotype","start_position","end_position",
                            "insideFeature","distancetoFeature","shortestDistance","fromOverlappingOrNearest")
diffbind_df_intersect_min[, (diffbind_cols_to_remove) := NULL]
# ~~~~~~~~~~~~~~~ #





# ~~~~~ CALCULATE LOG RATIO ~~~~~~~~ #
# create an empty data frame to fill with log ratio values
diff_log_ratio <- setNames(data.frame(matrix(nrow = nrow(diffbind_df_intersect_min), ncol = length(patientID_matches))),patientID_matches)
rownames(diff_log_ratio) <- diffbind_df_intersect_min[["external_gene_name"]]

# iterate over the patients; get only diff min cols per patient
for(i in seq_along(patientID_matches)){
  tmp_ID <- patientID_matches[i]
  # print(tmp_ID)
  # print(tmp_ID %in% colnames(diffbind_df_intersect_min))
  
  # get the columns with the patient DR
  tmp_colnames <- grep(pattern = tmp_ID,x = colnames(diffbind_df_intersect_min),value = TRUE)
  
  # if there aren't 2 colnames, break the loop and skip to next entry
  if (length(tmp_colnames) !=2) next
  
  # get the R colname
  tmp_colnames_R <- grep(pattern = paste0("^",tmp_ID,"[[:punct:]]R.*$"),x = tmp_colnames,perl = TRUE,value = TRUE)
  # get the D colname
  tmp_colnames_D <- grep(pattern = paste0("^",tmp_ID,"[[:punct:]]D.*$"),x = tmp_colnames,perl = TRUE,value = TRUE)
  
  # calculate the log ratio of the two columns
  diff_log_ratio[,tmp_ID] <- log2(diffbind_df_intersect_min[[tmp_colnames_R]] / diffbind_df_intersect_min[[tmp_colnames_D]])
  # diffbind_df_intersect_min[, .SD, .SDcols = tmp_colnames_R][[1]]
}

# remove GHW sample if present because its bad
if(any(grepl(pattern = "GHW",x = colnames(diff_log_ratio)))) diff_log_ratio <- diff_log_ratio[,grep(pattern = "GHW",x = colnames(diff_log_ratio),value = FALSE,invert = TRUE)]
# boxplot(diff_log_ratio)

# add a column in the data frame with the gene names, from the rownames
diff_log_ratio[["gene"]] <- rownames(diff_log_ratio)

# melt it into long format
library(reshape2)
diff_log_ratio_long <- reshape2::melt(diff_log_ratio,id.vars="gene",variable.name="sample",value.name="diff_peak_log_ratio")

# add a column to hold the 'status' e.g. Up regulated vs. Down regulated in the gene_expr table
diff_log_ratio_long[["gene_expression_status"]] <- NA
# ~~~~~~~~~~~~~~~ #




# ~~~~~ DETERMINE UP/DOWN REG GENE STATUS ~~~~~~~~ #
# add a column in the diff_log_ratio_long table telling whether each entry is up or down regulated in the gene expression dataset

# first make a subset of the LR columns
# # seach for the ones that match the diff_log_ratio samples
search_pattern2 <- paste(as.character(unique(diff_log_ratio_long[["sample"]])),collapse = "|")
# [1] "AGK|BVI|CBK ... 

# # search for matching colnames
colname_matches2 <- grep(pattern = search_pattern2,x = colnames(gene_expr_df_intersect),value = TRUE)

# # get the LR colnames to keep
LR_colnames <- grep(pattern = "^.*LR$",x = colname_matches2,perl = TRUE,value = TRUE)

# # make a df with just those columns
gene_expr_df_intersect_LR <- gene_expr_df_intersect[,LR_colnames]

# remove the .Exp.LR from the colnames..
colnames(gene_expr_df_intersect_LR) <- gsub(perl = TRUE,pattern = "^([[:alpha:]]*)([[:punct:]]*.*)$",replacement = "\\1",x = colnames(gene_expr_df_intersect_LR))

# add a column with the gene names for melting
gene_expr_df_intersect_LR[["gene"]] <- rownames(gene_expr_df_intersect_LR)


# melt it into long format
gene_expr_df_intersect_LR_long <- reshape2::melt(gene_expr_df_intersect_LR,id.vars="gene",variable.name="sample",value.name="gene_expression_log_ratio")
# ~~~~~~~~~~~~~~~ #


# ~~~~~ MERGE THE TABLES ~~~~~~~~ #
# merge together the two long format diff_log_ratio's and gene_expr log ratio's
diff_gene_expr_merged <- base::merge(gene_expr_df_intersect_LR_long,diff_log_ratio_long,by=c("gene","sample"))

# set the gene expression status to UP or DOWN based on gene expression value
# 1.5x up/down expression
diff_gene_expr_merged[["gene_expression_status"]] <- ifelse(diff_gene_expr_merged[["gene_expression_log_ratio"]]>=0.58, "UP",
                                                            ifelse(diff_gene_expr_merged[["gene_expression_log_ratio"]]<=-0.58,"DOWN",no = NA))
# ~~~~~~~~~~~~~~~ #


# # ~~~~~ write some files ~~~~~~~~ #
head(diff_gene_expr_merged)
# write.csv(file = paste0(main_outdir,"/diff_gene_expr_merged.csv"))
write.table(x = diff_gene_expr_merged,
            file = paste0(main_outdir,"/diff_gene_expr_merged.csv"),
            sep = '\t',quote = FALSE,col.names=NA,row.names = rownames(diff_gene_expr_merged))
# ~~~~~~~~~~~~~~~ #


# ~~~~~ MAKE PLOTS ~~~~~~~~ #
# make sure there are non-NA entries present; there's at least one UP or DOWN entry
if(nrow(diff_gene_expr_merged[which(diff_gene_expr_merged[["gene_expression_status"]]=="UP"),])>=50 &&
   nrow(diff_gene_expr_merged[which(diff_gene_expr_merged[["gene_expression_status"]]=="DOWN"),])>=50){
  # bean plot of all samples
  library(beanplot)
  pdf(file = paste0(main_outdir,"/all_samples.pdf"),height = 8,width = 12)
  beanplot(diff_peak_log_ratio~gene_expression_status*sample,
           data=diff_gene_expr_merged,
           what=c(0,1,1,0),
           border = NA,
           col=list('grey','blue'),
           ylab="DiffBind ratio = log2 ( Diff peak R / Diff peak D )",
           main=paste0("DiffBind results for ",histone_mark," ",params_branch), cex.main=0.9
           ,side = "both"
  )
  legend('topright', fill=c('grey','blue'), legend= c('Down Genes', 'Up Genes'))
  abline(h=0,col="red")
  dev.off()
  
}
# ~~~~~~~~~~~~~~~ #






