#!/usr/bin/env Rscript

## USAGE: multi_peaks_UpSet_plot.R <sampleID> /path/to/venn.txt
## This script will process venn output from HOMER mergePeaks and make an UpSet plot
## this is currently only for R version 3.3.0; module unload r; module load r/3.3.0; R

# resources:
# http://www.caleydo.org/tools/upset/
# https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html#example-5-alternative-input-formats 
# https://github.com/hms-dbmi/UpSetR
 
# ~~~~~ LOAD PACKAGES ~~~~~~~ #
library("UpSetR")
library("ggplot2")
library("grid")
library("plyr")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 


# ~~~~~ GET SCRIPT ARGS ~~~~~~~ #
args <- commandArgs(TRUE)
cat("Script args are:\n")
args

SampleID<-args[1]
venn_table_file<-args[2]

plot_outdir<-dirname(venn_table_file)
plot_filepath<-file.path(plot_outdir, paste0(SampleID,"_UpSetR_plot.pdf"))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 


# ~~~~~ PARSE THE VENN TABLE ~~~~~~~ #
# read in the venn text
venn_table_df <- read.table(venn_table_file,header = TRUE,sep = "\t",stringsAsFactors = FALSE,check.names = FALSE)
# > venn_table_df
# AGK-D-H3K27AC.bed AGK-R-H3K27AC.bed Total                                Name
# 1                                   X 29005                   AGK-R-H3K27AC.bed
# 2                 X                   19523                   AGK-D-H3K27AC.bed
# 3                 X                 X 35344 AGK-D-H3K27AC.bed|AGK-R-H3K27AC.bed


# get the venn categories
venn_categories<-colnames(venn_table_df)[!colnames(venn_table_df) %in% c("Total","Name")] 
cat("Venn categories are:\n")
venn_categories
# > cat("Venn categories are:\n"); venn_categories
# Venn categories are:
#   [1] "AGK-D-H3K27AC.bed" "AGK-R-H3K27AC.bed"


# venn_categories
num_categories <- length(venn_categories)
cat("Num categories are:\n")
num_categories
# > num_categories<-length(venn_categories)
# > cat("Num categories are:\n"); num_categories
# Num categories are:
#   [1] 2


# make a summary table
venn_summary<-venn_table_df[!colnames(venn_table_df) %in% venn_categories]
cat("Venn summary table is categories are:\n")
venn_summary
# > venn_summary<-venn_table_df[!colnames(venn_table_df) %in% venn_categories]
# > cat("Venn summary table is categories are:\n"); venn_summary
# Venn summary table is categories are:
#   Total                                Name
# 1 29005                   AGK-R-H3K27AC.bed
# 2 19523                   AGK-D-H3K27AC.bed
# 3 35344 AGK-D-H3K27AC.bed|AGK-R-H3K27AC.bed
 

# write summary table
# plot_filepath<-paste0(plot_outdir,"/",SampleID,"_UpSetR_plot.pdf") 
# write.table(venn_summary,file = paste0(plot_outdir,"/","venn_summary.tsv"),quote = FALSE,row.names = FALSE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 


# ~~~~~ SET UP THE PLOT ~~~~~~~ #
# convert the summary table to a numeric/int vector, with element names as the combination names
# # swap the | character with &; for passing to UpSet fromExpression
upset_expression <- setNames(venn_summary[['Total']], gsub("|","&",venn_summary[['Name']], fixed = TRUE))
# > upset_expression
# AGK-R-H3K27AC.bed                   AGK-D-H3K27AC.bed
# 29005                               19523
# AGK-D-H3K27AC.bed&AGK-R-H3K27AC.bed
# 35344

# save the plot into a PDF
cat("Plot output file is:\n")
plot_filepath

pdf(plot_filepath, width = 8, height = 8, onefile=FALSE)
upset(fromExpression(upset_expression), 
      nsets = num_categories, 
      order.by = "freq", 
      decreasing = T, 
      mainbar.y.label = "Overlapping Peaks", 
      sets.x.label = "Peaks per Category") # , group.by = "sets"
grid.text(SampleID, vp = viewport(layout.pos.row = 2, layout.pos.col = 1), just = c("left", "top"))
dev.off()

