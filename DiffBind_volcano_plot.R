#!/usr/bin/env Rscript

# this script contains functions to parse DiffBind ChIP-Seq differential binding data and create volcano plots

# function to create volcano plot from DiffBind output CSV file
DiffBind_volcano_plot <- function(diffbind_file, signif_p_value=0.05, fold_change_value=1, plot_colors=c("black","red","orange","blue")){
	
	# read in file
	diff_df <- read.delim(file = diffbind_file,header = TRUE,sep = ',')
	
	# get sample ID's
	sampleID_1 <- colnames(diff_df)[7]
	sampleID_2 <- colnames(diff_df)[8]
	
	# signif_p_value <- 0.05
	# fold_change_value <- 1
	
	# Make a basic volcano plot
	# plot_colors <- c("black","red","orange","blue")
	
	# set layout matrix
	plot_layout_matrix<-structure(c(1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 1L, 2L, 
																	2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 
																	3L, 3L, 3L, 3L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L),
																.Dim = c(11L,4L), 
																.Dimnames = list(NULL, c("V1", "V2", "V3", "V4")))
	
	# subset the matrix for two panels # keep this because I like this code sample
	plot_layout_matrix2 <- plot_layout_matrix[(rowSums(plot_layout_matrix < 3) > 0), , drop = FALSE]
	# x[which(x < 3, arr.ind = TRUE)[,1],]
	
	
	# setup the panel layout
	layout(plot_layout_matrix2) 
	# adjust first panel margins
	# plot margins # c(bottom, left, top, right) # default is c(5, 4, 4, 2) + 0.1
	par(mar=c(0,6,7,0))
	# call blank plot to fill the first panel
	plot(1,type='n',axes=FALSE,xlab="",ylab="",main = paste0("Volcano plot: ", sampleID_1, " vs. ",sampleID_2),cex.main=1) 
	# set up the Legend in the first panel
	legend("bottom",legend=c("Insignificant",
													 paste0("p < ",signif_p_value),
													 paste0("Fold Change > ",fold_change_value),
													 paste0("p < ",signif_p_value," & Fold Change > ",fold_change_value)),
				 fill=plot_colors,bty = "n",ncol=length(plot_colors),cex=0.9)
	# adjust second panel margins
	par(mar=c(6,4,0,3)+ 0.1)
	# start second panel plots
	with(diff_df, plot(Fold, -log10(p.value), pch=20, xlim=c(min(Fold)-1,max(Fold)+1),col=plot_colors[1],xlab = "Fold Change"))
	# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
	with(subset(diff_df, p.value<signif_p_value ), points(Fold, -log10(p.value), pch=20, col=plot_colors[2]))
	with(subset(diff_df, abs(Fold)>fold_change_value), points(Fold, -log10(p.value), pch=20, col=plot_colors[3]))
	with(subset(diff_df, p.value<signif_p_value & abs(Fold)>fold_change_value), points(Fold, -log10(p.value), pch=20, col=plot_colors[4]))
	
	# Label points with the textxy function from the calibrate plot
	# install.packages("calibrate")
	# library("calibrate")
	# with(subset(diff_df, p.value<signif_p_value & abs(Fold)>fold_change_value), textxy(Fold, -log10(p.value), labs=external_gene_name, cex=.8))
	
}

# ~~~~~~~~~~~~ # ~~~~~~~~~~~~~~ # 
out_dir <- "/ifs/home/kellys04/projects/SmithLab_ChIPSeq-RNASeq_2016-12-31/project_notes/diffbind_files"

diff_cl4_ctrl_file <- paste0(out_dir,"/diff_bind.Sample4-ChIPSeq-vs-Control-ChIPSeq.p100.csv")
diff_cl4_cl5_file <- paste0(out_dir,"/diff_bind.Sample4-ChIPSeq-vs-Sample5-ChIPSeq.p100.csv")
diff_cl5_ctrl_file <- paste0(out_dir,"/diff_bind.Sample5-ChIPSeq-vs-Control-ChIPSeq.p100.csv")

pdf(file = paste0(out_dir,"/volcano_plot_Sample4-ChIPSeq-vs-Control-ChIPSeq.p100.pdf"),height = 8,width = 8)
DiffBind_volcano_plot(diffbind_file = diff_cl4_ctrl_file)
dev.off()

pdf(file = paste0(out_dir,"/volcano_plot_Sample4-ChIPSeq-vs-Sample5-ChIPSeq.p100.pdf"),height = 8,width = 8)
DiffBind_volcano_plot(diffbind_file = diff_cl4_cl5_file)
dev.off()

pdf(file = paste0(out_dir,"/volcano_plot_Sample5-ChIPSeq-vs-Control-ChIPSeq.p100.pdf"),height = 8,width = 8)
DiffBind_volcano_plot(diffbind_file = diff_cl5_ctrl_file)
dev.off()
# ~~~~~~~~~~~~ # ~~~~~~~~~~~~~~ # 


DiffBind_volcano_plot_top_promoters <- function(diffbind_file, signif_p_value=0.05, fold_change_value=1, plot_colors=c("grey","red","orange","blue", "green"), top_gene_number=10){
	# this function creates a volcano plot directly from a DiffBind .csv file
	# only plots promoter genes from DiffBind output
	# highlights top x fold changed genes
	
	# ~~~~~~ PROCESS DIFFBIND FILE ~~~~~~~~ # 
	# read in file
	diff_df <- read.delim(file = diffbind_file,header = TRUE,sep = ',')

	
	# get the colname order for later
	original_colnames <- colnames(diff_df)
	
	library("data.table")
	# GET THE CLOSEST DIFFBIND PEAK PER GENE 
	# in Diffbind subset, getthe peaks closest to each gene
	# # for each gene name, find the peak with the lowest "shortestDistance" value
	
	# make it a data.table type
	setDT(diff_df)
	# To get the promoter genes, use only the DiffBind peaks with the shortest distance per gene;
	# get min value per factor, also remove duplicates
	diff_df_min <- as.data.frame(diff_df[, .SD[which.min(shortestDistance)], by=external_gene_name])
	# fix the colname order
	diff_df_min <- diff_df_min[,original_colnames]
	
	# get sample ID's
	sampleID_1 <- colnames(diff_df_min)[7]
	sampleID_2 <- colnames(diff_df_min)[8]
	
	
	# subset for significant genes
	diff_df_min_signiff <- as.data.frame(subset(diff_df_min, p.value<signif_p_value & abs(Fold)>fold_change_value))
	
	
	# subset top up/down genes
	diff_df_min_signiff_UpFold <- diff_df_min_signiff[with(diff_df_min_signiff, order(Fold, p.value)), c("Fold","p.value","external_gene_name")][1:top_gene_number,]
	
	diff_df_min_signiff_DnFold <- diff_df_min_signiff[with(diff_df_min_signiff, order(-Fold, p.value)), c("Fold","p.value","external_gene_name")][1:top_gene_number,]
	# ~~~~~~~~~~~~~~ # 
	
	
	# ~~~~~~ START PLOT ~~~~~~~~ # 
	# set the multi-pane plot layout matrix
	plot_layout_matrix<-structure(c(1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 1L, 2L, 
																	2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 
																	3L, 3L, 3L, 3L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L),
																.Dim = c(11L,4L), 
																.Dimnames = list(NULL, c("V1", "V2", "V3", "V4")))
	
	# subset the matrix for two panels # keep this because I like this code sample
	plot_layout_matrix2 <- plot_layout_matrix[(rowSums(plot_layout_matrix < 3) > 0), , drop = FALSE]
	# x[which(x < 3, arr.ind = TRUE)[,1],]
	
	# START PLOT
	# setup the panel layout
	layout(plot_layout_matrix2) 
	
	# FIRST PLOT PANEL
	# adjust first panel margins
	# plot margins # c(bottom, left, top, right) # default is c(5, 4, 4, 2) + 0.1
	par(mar=c(1,3,1,0))
	# call blank plot to fill the first panel
	plot(1,type='n',axes=FALSE,xlab="",ylab="")
	# set up the Legend in the first panel; 1st legend only title, 2nd legend only legend
	legend("center",legend = "",title=paste0("Volcano plot: ", sampleID_1, " vs. ",sampleID_2,"\n"),cex=1.3, bty='n') 
	legend("bottom",legend=c("Insignificant",
													 paste0("p < ",signif_p_value),
													 paste0("log2 Fold Change > ",fold_change_value),
													 paste0("p < ",signif_p_value," & Fold Change > ",fold_change_value),
													 paste0("Top log2 Fold Change genes")),
				 fill=plot_colors,bty = "n",ncol=2,cex=0.9)
	
	# SECOND PLOT PANEL
	# adjust second panel margins
	par(mar=c(6,4,0,3)+ 0.1)
	
	# add base volcano plot
	with(diff_df_min, plot(Fold, -log10(p.value), pch=20, xlim=c(min(Fold)-1,max(Fold)+1),col=plot_colors[1],xlab = "Fold Change"))
	
	# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
	with(subset(diff_df_min, p.value<signif_p_value ), points(Fold, -log10(p.value), pch=20, col=plot_colors[2]))
	with(subset(diff_df_min, abs(Fold)>fold_change_value), points(Fold, -log10(p.value), pch=20, col=plot_colors[3]))
	with(subset(diff_df_min, p.value<signif_p_value & abs(Fold)>fold_change_value), points(Fold, -log10(p.value), pch=20, col=plot_colors[4]))
	
	# add points and gene labels for top genes
	library("calibrate")
	with(na.omit(diff_df_min_signiff_UpFold), points(Fold, -log10(p.value), pch=20, col="green"))
	with(na.omit(diff_df_min_signiff_UpFold), textxy(Fold, -log10(p.value), labs=external_gene_name, cex=.8))

	with(na.omit(diff_df_min_signiff_DnFold), points(Fold, -log10(p.value), pch=20, col="green"))
	with(na.omit(diff_df_min_signiff_DnFold), textxy(Fold, -log10(p.value), labs=external_gene_name, cex=.8))
	
}
	



