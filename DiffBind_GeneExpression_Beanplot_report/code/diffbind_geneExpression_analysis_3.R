#!/usr/bin/env Rscript

# !! third attempt at this pipeline to try and make things clearer from the start !!
# this time make some nice ggplot2 split density plots and violin plots along with beanplots

## USAGE: diffbind_geneExpression_analysis.R /path/to/outdir /path/to/microarray_gene_expression.tsv /path/to/diffbind.csv <HistMark> <analysis_params_branch>

## OUTPUT: <outdir>/<SampleID>/<HistMark>/<analysis_params_branch>/[paired_boxplot.pdf genes_up.tsv genes_down.tsv]
## e.g. <outdir>/ZNK/INPUT/peaks.by_sample.macs_broad-cutoff0025

## DESCRIPTION:
# This pipeline script will compare microarray gene expression data with
# DiffBind differential peak binding data
# the script will find the set of genes that are up or down regulated in the microarray data
# and create subsets of genes with differentially bound peaks from the DiffBind data
# and plot the DiffBind values for both subsets (up and down regulated genes) in a set of beanplots
# the subset diffbind values will be output as TSV, and the plot

library("data.table")
library("beanplot")
library("ggplot2")
library(reshape2)

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
# "H3K27AC"

# get the analysis params branch for the samples
params_branch <- args[5]
# "peaks.by_sample.macs_broad"
# ~~~~~~~~~~~~~~~ #



# ~~~~~~ READ IN FILE ~~~~~~~~ #
# diffbind data
diffbind_df <- read.csv(diffbind_file)
# gene expression data
gene_expr_df <- read.table(gene_expr_file,sep = '\t',header = TRUE,row.names = 1,quote = "")
# ~~~~~~~~~~~~~~~ #


# ~~~~~ GET SAMPLE IDS ~~~~~~~~ #
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





# ~~~~~ CALCULATE LOG RATIO FULL DIFFBIND DATASET ~~~~~~~~ #
# create an empty data frame to fill with log ratio values
diff_log_ratio_full <- setNames(data.frame(matrix(nrow = nrow(diffbind_df), ncol = length(patientID_matches))),patientID_matches)
# add the genes from DiffBind
diff_log_ratio_full[["gene"]] <- diffbind_df[["external_gene_name"]]

# enter the log ratio values into the df per patient
# # iterate over the patients; get only diff min cols per patient
for(i in seq_along(patientID_matches)){
  tmp_ID <- patientID_matches[i]
  print(tmp_ID)
  # print(tmp_ID %in% colnames(diffbind_df_intersect_min))
  
  # get the diffBind columns with the patient's ID
  tmp_colnames <- grep(pattern = tmp_ID,x = colnames(diffbind_df),value = TRUE)
  
  # make sure there are 2 colnames, otherwise skip to the next iteration of the loop!
  if (length(tmp_colnames) !=2) next
  
  # get the R colname
  tmp_colnames_R <- grep(pattern = paste0("^",tmp_ID,"[[:punct:]]R.*$"),x = tmp_colnames,perl = TRUE,value = TRUE)
  
  # get the D colname
  tmp_colnames_D <- grep(pattern = paste0("^",tmp_ID,"[[:punct:]]D.*$"),x = tmp_colnames,perl = TRUE,value = TRUE)
  
  # calculate the log ratio of the two columns
  # # log2 ( R / D )
  diff_log_ratio_full[,tmp_ID] <- log2(diffbind_df[[tmp_colnames_R]] / diffbind_df[[tmp_colnames_D]])
  
}

# remove GHW sample if present because its bad
if(any(grepl(pattern = "GHW",x = colnames(diff_log_ratio_full)))) diff_log_ratio_full <- diff_log_ratio_full[,grep(pattern = "GHW",x = colnames(diff_log_ratio_full),value = FALSE,invert = TRUE)]
# boxplot(diff_log_ratio_full)

# melt it into long format
diff_log_ratio_full_long <- reshape2::melt(diff_log_ratio_full,id.vars="gene",variable.name="sample",value.name="diff_peak_log_ratio")

# add a column to hold the 'status' e.g. Up regulated vs. Down regulated in the gene_expr table
diff_log_ratio_full_long[["gene_expression_status"]] <- NA
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

# make it a data.table type
setDT(diffbind_df_intersect)
# get min value per factor, also remove duplicates
diffbind_df_intersect_min <- diffbind_df_intersect[, .SD[which.min(shortestDistance)], by=external_gene_name]

# remove extraneous columns
diffbind_cols_to_remove <-c("seqnames","start","end","width","strand","Conc","Conc_D","Conc_R","Fold",
                            "p.value","FDR","feature","gene_biotype","start_position","end_position",
                            "insideFeature","distancetoFeature","shortestDistance","fromOverlappingOrNearest")
diffbind_df_intersect_min[, (diffbind_cols_to_remove) := NULL]
# ~~~~~~~~~~~~~~~ #





# ~~~~~ DETERMINE UP/DOWN REG GENE STATUS ~~~~~~~~ #
# add a column in the diff_log_ratio_long table telling whether each entry is up or down regulated in the gene expression dataset

# first make a subset of the LR columns
# # seach for the ones that match the diff_log_ratio samples
search_pattern2 <- paste(as.character(unique(diff_log_ratio_full_long[["sample"]])),collapse = "|")
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
diff_gene_expr_merged <- base::merge(gene_expr_df_intersect_LR_long,diff_log_ratio_full_long,by=c("gene","sample"),all=TRUE)

# set the gene expression status to UP or DOWN based on gene expression value
# 1.5x up/down expression
diff_gene_expr_merged[["gene_expression_status"]] <- ifelse(diff_gene_expr_merged[["gene_expression_log_ratio"]]>=0.58, "UP",
                                                            ifelse(diff_gene_expr_merged[["gene_expression_log_ratio"]]<=-0.58,"DOWN",no = NA))
# ~~~~~~~~~~~~~~~ #
# 

# # ~~~~~ write some files ~~~~~~~~ #
head(diff_gene_expr_merged)
# write.csv(file = paste0(main_outdir,"/diff_gene_expr_merged.csv"))
write.table(x = diff_gene_expr_merged,
            file = paste0(main_outdir,"/diff_gene_expr_merged.tsv"),
            sep = '\t',quote = FALSE,col.names=NA,row.names = rownames(diff_gene_expr_merged))
# ~~~~~~~~~~~~~~~ #





# ~~~~~ PLOT DIFFBIND VALUES ~~~~~~~~ #
# ggplot2 basic beanplot
theme_set(theme_bw()) # use the black and white theme throughout
# base violin plot
p <- ggplot(diff_gene_expr_merged, aes(x=sample, y=diff_peak_log_ratio, group=sample)) + geom_violin(trim=FALSE,fill="darkorchid4",color=NA,scale="width")
# change the color theme
p <- p + theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank()) 
# set the axis labels
p <- p + xlab("") + ylab("DiffBind ratio = log2 ( Diff peak R / Diff peak D )") 
# add chart title
p <- p + ggtitle(paste0("DiffBind results for ",histone_mark," ",params_branch))


# overlay with points based on the gene expression data
# # hold off on these, the jitter plot especially does not come out well but this code might be useful
# k <- ggplot(diff_gene_expr_merged, aes(x=sample, y=diff_peak_log_ratio, group=sample)) + geom_violin(trim=FALSE,fill="darkorchid4",color=NA)
# q <- p + geom_point(data=na.omit(diff_gene_expr_merged), aes(colour = gene_expression_status),position = position_jitter(width = 0.3, height = 0))
# r <- p + geom_jitter(data=na.omit(diff_gene_expr_merged), alpha=I(1/4), aes(color=gene_expression_status),width = 0.1, height = 0)
# print(q)
# print(p) # base beanplot, full width
# print(k) # base beanplot, width = count
# print(r) # full width beanplot with jittered points overlayed for gene expression

# ggplot2 split density plots; these actually look nicer than the beanplots and are more flexible in ggplot2
# start the density plot
dens_adjust=1 # plot denisty adjustment value; default = 1
tt <- ggplot(data=diff_gene_expr_merged)
# add up genes
tt <- tt + geom_density(data=subset(diff_gene_expr_merged,gene_expression_status == unique(diff_gene_expr_merged[["gene_expression_status"]])[2]), aes(x=diff_peak_log_ratio,y=..density..,fill=gene_expression_status), trim=F,alpha=0.4,adjust = dens_adjust,color=NA)
# add down genes
tt <- tt + geom_density(data=subset(diff_gene_expr_merged,gene_expression_status == unique(diff_gene_expr_merged[["gene_expression_status"]])[3]), aes(x=diff_peak_log_ratio,y=..density..,fill=gene_expression_status), trim=F,alpha=0.4,adjust = dens_adjust,color=NA)
# add all genes
tt <- tt + geom_density(data=diff_gene_expr_merged, aes(x=diff_peak_log_ratio,y=-..density..,fill="ALL", trim=F),alpha=1,adjust = dens_adjust,color=NA)
# rotate axis
tt <- tt + coord_flip() 
# add x labels
tt <- tt + xlab("DiffBind ratio = log2 ( Diff peak R / Diff peak D )") 
# set the plot faceting to split by samples
tt <- tt + facet_grid(~sample) 
# set the plot title
tt <- tt + ggtitle(paste0("DiffBind results for ",histone_mark," ",params_branch))
# set the theme attributes
tt <- tt + theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x =  element_blank(), legend.box = "horizontal", legend.position = "bottom", strip.text = element_text(face = "bold"), plot.title = element_text(margin=margin(b = 10,t = 10, unit = "pt"),face = "bold"), strip.background = element_rect(fill="lightgrey")) 
# fix the legend title
tt <- tt + labs(fill="Genes: ") 
tt <- tt + scale_fill_manual( values = c("darkorchid4","blue","red"))
# print(tt)

# Only make a plot if there are at least 50 genes in both UP and DOWN statuses; 
# # I might change these labels later so don't hard code them!
# > unique(diff_gene_expr_merged[["gene_expression_status"]])
# [1] NA     "DOWN" "UP"  
# if you don't run this test you get a bunch of wacky looking bad plots and some empty canvas PDFs which break the report
if(nrow(diff_gene_expr_merged[which(diff_gene_expr_merged[["gene_expression_status"]]==unique(diff_gene_expr_merged[["gene_expression_status"]])[2]),])>=50 &&
   nrow(diff_gene_expr_merged[which(diff_gene_expr_merged[["gene_expression_status"]]==unique(diff_gene_expr_merged[["gene_expression_status"]])[3]),])>=50){
  # base R beanplot
  pdf(file = paste0(main_outdir,"/beanplot_all_samples1.pdf"),height = 8,width = 12)
  beanplot(diff_peak_log_ratio~sample,
           # data=diff_log_ratio_full_long,
           data=diff_gene_expr_merged,
           what=c(0,1,1,0),
           # what=c(0,1,1,1),
           border = NA,
           bw="nrd",
           # col=list('grey','purple'),
           col='purple',
           ylab="DiffBind ratio = log2 ( Diff peak R / Diff peak D )",
           main=paste0("DiffBind results for ",histone_mark," ",params_branch), cex.main=0.9
           # ,side = "both"
  )
  # hold off on this strip chart overlay:
  # stripchart(diff_peak_log_ratio~gene_expression_status*sample,data = diff_gene_expr_merged,
  #            vertical=TRUE,add=TRUE,cex=0.1,method="jitter",jitter=0.05)
  dev.off()
  
  # ggplot2 violin plot
  pdf(file = paste0(main_outdir,"/ggbeanplot_plain.pdf"),height = 8,width = 12)
  print(p)
  dev.off()
  
  # ggplot2 split density plot
  pdf(file = paste0(main_outdir,"/ggbeanplot_split.pdf"),height = 8,width = 12)
  print(tt)
  dev.off()
}

# ~~~~~~~~~~~~~~~ #
