#!/usr/bin/env Rscript

# a quick short pipeline to run Gene Ontology analysis on the DiffBind differential peak binding output
# make some plots showing the GO results

# tested with Bioconductor 3.3 & clusterProfiler 3.0
# source("https://bioconductor.org/biocLite.R")
# biocLite("clusterProfiler")
# biocLite("org.Hs.eg.db")
# biocLite("AnnotationHub")
library("clusterProfiler") # https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html
library("org.Hs.eg.db") # keytypes(org.Hs.eg.db)
library("AnnotationHub")



out_dir <- "/ifs/home/kellys04/projects/SmithLab_ChIPSeq-RNASeq_2016-12-31/project_notes/gene_ontology"
diff_input_dir <- "/ifs/home/kellys04/projects/SmithLab_ChIPSeq-RNASeq_2016-12-31/project_notes/diffbind_files"

diff_cl4_ctrl_file <- paste0(diff_input_dir,"/diff_bind.Sample4-vs-Control-ChIPSeq.p100.csv")



# ~~~~~ GET DIFFBIND GENES~~~~~~~~~ #
# read in diffbind file
diff_df <- read.delim(file = diff_cl4_ctrl_file,header = TRUE,sep = ',')

# subset for significant entries
signif_p_value <- 0.05
diff_df <- subset(diff_df, p.value<signif_p_value)


# peel off the Fold change values, named by gene name
diff_fold_genes <- diff_df[,c("Fold","external_gene_name")]
# > colnames(diff_fold_genes)
# [1] "Fold"               "external_gene_name"


# Get EntrezID's for the genes; returns a new data frame
diff_entrezID <- bitr(diff_fold_genes[["external_gene_name"]], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# > colnames(diff_entrezID)
# [1] "SYMBOL"   "ENTREZID"

# merge the sets
# # drops genes that didn't have EntrezID's, e.g. snoRNA, misc RNA's, non-coding, etc
diff_fold_genes_ID <- merge(diff_fold_genes,diff_entrezID,by.x="external_gene_name", by.y="SYMBOL") # ,all=TRUE
# > colnames(diff_fold_genes_ID)
# [1] "external_gene_name" "Fold"               "ENTREZID"

# peel off a vector of fold changes, named by EntrezID
fold_change <- setNames(diff_fold_genes_ID[["Fold"]],diff_fold_genes_ID[["ENTREZID"]])

# drop fold change entries less than 2
fold_change_cutoff <- 2
fold_change <- fold_change[fold_change >= fold_change_cutoff]




# ~~~~~ GO ANALYSIS ~~~~~~~~~ #
# from the vignette
# 5.2 GO classification
ggo <- groupGO(gene     = names(fold_change),
							 OrgDb    = org.Hs.eg.db,
							 # organism="human",
							 ont      = "BP",
							 level    = 3,
							 readable = TRUE)
head(summary(ggo))

# 5.3 GO over-representation test
ego <- enrichGO(gene          = names(fold_change), # genes to look for enrichment in
								# universe      = names(GO_genes), # background genes
								OrgDb         = org.Hs.eg.db,
								# organism="human",
								ont           = "BP",
								pAdjustMethod = "BH",
								pvalueCutoff  = 0.01,
								qvalueCutoff  = 0.05, 
								readable      = TRUE)
head(summary(ego))


# plots
pdf(file = paste0(out_dir,"/GO_plots.pdf"),width = 10,height = 10,onefile = TRUE)
barplot(ggo, drop=TRUE, showCategory=12)
barplot(ego, showCategory=8)
dotplot(ego)
enrichMap(ego,vertex.label.cex=0.7)
cnetplot(ego, categorySize="pvalue", foldChange=fold_change,vertex.label.cex=0.7)
cnetplot(ego, categorySize="geneNum", foldChange=fold_change,vertex.label.cex=0.7)
dev.off()



# ~~~~~ GO CUSTOM WORKFLOW ~~~~~~~~~ #
# custom GO functions
custom_groupGO <- function(gene_list, sub_ontology="BP"){
	# gene_list: char vector of gene ID's
	# sub_ontology # One of "MF", "BP", and "CC" subontologies.; Cellular Component, ... ... 
	
	# run the group GO analysis
	ggo <- groupGO(gene     = gene_list,
								 OrgDb    = org.Hs.eg.db,
								 ont      = sub_ontology,
								 level    = 3,
								 readable = TRUE)
	return(ggo)
}


custom_enrichGO <- function(gene_list, sub_ontology="BP"){
	# gene_list: char vector of gene ID's
	# sub_ontology # One of "MF", "BP", and "CC" subontologies.
	
	# run the enrichment GO analysis
	ego <- enrichGO(gene          = gene_list, # genes to look for enrichment in
									# universe      = names(GO_genes), # background genes
									OrgDb         = org.Hs.eg.db,
									ont           = sub_ontology,
									pAdjustMethod = "BH",
									pvalueCutoff  = 0.01,
									qvalueCutoff  = 0.05, 
									readable      = TRUE)
	return(ego)
}



get_GO_results <- function(GO_object){
	# save the results to a new dataframe
	ggo_results <- GO_object@result
	
	# sort by count
	ggo_results <- ggo_results[with(ggo_results, order(Count,decreasing = TRUE)),] # 
	return(ggo_results)
}


custom_groupGO_barplot <- function(GO_counts, GO_descriptions, sub_ontology="BP"){
	# make a vector of values to plot
	GO_results_vec <- setNames(GO_counts,GO_descriptions)
	
	plot_cols <- rainbow(n = length(GO_results_vec))
	
	# plot margins # c(bottom, left, top, right) # default is c(5, 4, 4, 2) + 0.1
	par(mar=c(5, 20, 4, 2) + 0.1)
	barplot(GO_results_vec,horiz = TRUE,las=1,col = plot_cols,xlab = "Number of genes",cex.names=0.75,main = paste0("Group GO plot: ",sub_ontology))
}


custom_enrichGO_barplot <- function(GO_counts, GO_descriptions, sub_ontology="BP"){
	# make a vector of values to plot
	GO_results_vec <- setNames(GO_counts,GO_descriptions)
	
	plot_cols <- rainbow(n = length(GO_results_vec))
	
	# plot margins # c(bottom, left, top, right) # default is c(5, 4, 4, 2) + 0.1
	par(mar=c(5, 23, 4, 2) + 0.1)
	barplot(GO_results_vec,horiz = TRUE,las=1,col = plot_cols,xlab = "Number of genes",cex.names=0.75,main = paste0("Enrichment GO plot: ",sub_ontology))
}



# save custom GO barplots to a single PDF
pdf(file = paste0(out_dir,"/GO_plots5.pdf"),width = 10,height = 10,onefile = TRUE)
for(sub_ont in c("MF","BP","CC")){
	# get the subontology category
	tmp_ont <- sub_ont
	
	# run group GO functions
	tmp_ggo <- custom_groupGO(names(fold_change), tmp_ont)
	tmp_ggo_results <- get_GO_results(tmp_ggo)
	# save results to table
	write.table(x = tmp_ggo_results, file = paste0(out_dir,"/",tmp_ont,"_ggo.tsv"),quote = FALSE,sep = '\t',row.names = FALSE,col.names = TRUE)
	# plot
	custom_groupGO_barplot(tmp_ggo_results[["Count"]][1:10],tmp_ggo_results[["Description"]][1:10],tmp_ont)
	
	# run enrichment GO functions
	tmp_ego <- custom_enrichGO(names(fold_change), tmp_ont)
	tmp_ego_results <- get_GO_results(tmp_ego)
	# save results to table
	write.table(x = tmp_ego_results, file = paste0(out_dir,"/",tmp_ont,"_ego.tsv"),quote = FALSE,sep = '\t',row.names = FALSE,col.names = TRUE)
	# plot
	custom_enrichGO_barplot(tmp_ego_results[["Count"]][1:10],tmp_ego_results[["Description"]][1:10],tmp_ont)
}
dev.off()
