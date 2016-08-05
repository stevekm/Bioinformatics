#!/usr/bin/env Rscript

# a quick short pipeline to run Gene Ontology analysis on the DiffBind differential peak binding output
# make some plots showing the GO results
# this includes the standard analysis steps as per the vignette, 
# and also custom functions I wrote to make the process easier

# tested with Bioconductor 3.3 & clusterProfiler 3.0
# source("https://bioconductor.org/biocLite.R")
# biocLite("clusterProfiler")
# biocLite("org.Hs.eg.db")
# biocLite("AnnotationHub")
library("clusterProfiler") # https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html
library("org.Hs.eg.db") # keytypes(org.Hs.eg.db)
library("AnnotationHub")



# ~~~~~ STANDARD WORKFLOW ~~~~~~~~~ # 
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
# https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html

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


# 6.1 KEGG over-representation test
kk <- enrichKEGG(gene         = names(fold_change), # vector of gene Entrez ID's
								 organism     = 'hsa',
								 pvalueCutoff = 0.05)
head(summary(kk))
cnetplot(kk, categorySize="geneNum", foldChange=geneList)

# 6.2 KEGG Gene Set Enrichment Analysis

kk2 <- gseKEGG(geneList     = fold_change, # vector of Fold Changes named by gene ID
							 organism     = 'hsa',
							 nPerm        = 1000,
							 minGSSize    = 120,
							 pvalueCutoff = 0.05,
							 verbose      = FALSE)
head(summary(kk2))
gseaplot(kk2, geneSetID = "hsa04145")


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
# ~~~~~ CUSTOM FUNCTIONS ~~~~~~~~~ # 

# process the DiffBind data
get_diffbind_gene_foldchange <- function(diffbind_df, signif_p_value=0.05, diff_fold_labl="Fold", diff_gene_label="external_gene_name", fold_change_cutoff=2, sort_decreasing=FALSE){
	# diffbind_df : data frame containing diffbind data
	# signif_p_value : numeric pvalue cutoff, or FALSE for no cutoff
	# fold_change_cutoff : numeric fold change cutoff, or FALSE
	# diff_fold_labl : character label of the fold change column name in the diffbind dataframe
	# diff_gene_label : character label of the gene column name in the diffbind dataframe
	# sort_decreasing : logical, should the fold change values be sorted? 
	
	# read in diffbind file
	# diff_df <- read.delim(file = diff_cl4_ctrl_file,header = TRUE,sep = ',')
	
	# subset for significant entries
	if(signif_p_value) diffbind_df <- subset(diffbind_df, p.value<signif_p_value)
	
	
	
	
	# peel off the Fold change values, named by gene name
	diff_fold_genes <- diffbind_df[,c(diff_fold_labl,diff_gene_label)]
	# > colnames(diff_fold_genes)
	# [1] "Fold"               "external_gene_name"
	
	
	# Get EntrezID's for the genes; returns a new data frame
	diff_entrezID <- bitr(diff_fold_genes[[diff_gene_label]], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
	# > colnames(diff_entrezID)
	# [1] "SYMBOL"   "ENTREZID"
	
	# merge the sets
	# # drops genes that didn't have EntrezID's
	diff_fold_genes_ID <- merge(diff_fold_genes,diff_entrezID,by.x=diff_gene_label, by.y="SYMBOL") # ,all=TRUE
	# > colnames(diff_fold_genes_ID)
	# [1] "external_gene_name" "Fold"               "ENTREZID"
	
	# peel off a vector of fold changes, named by EntrezID
	fold_change <- setNames(diff_fold_genes_ID[[diff_fold_labl]],diff_fold_genes_ID[["ENTREZID"]])
	
	# drop fold change entries less than 2
	if(fold_change_cutoff) fold_change <- fold_change[fold_change >= fold_change_cutoff]
	
	# sort the fold change values
	# fold_change[order(fold_change,decreasing = TRUE)]
	if(sort_decreasing) fold_change <- fold_change[order(fold_change,decreasing = TRUE)]
	
	return(fold_change)
}

# run the GO tests
custom_groupGO <- function(gene_list, sub_ontology="BP", results_only=FALSE, sort_results=FALSE, sort_results_by="Count"){
	# gene_list: char vector of gene ID's
	# sub_ontology # One of "MF", "BP", and "CC" subontologies.; Cellular Component, ... ... 
	
	# run the group GO analysis
	ggo <- groupGO(gene     = gene_list,
								 OrgDb    = org.Hs.eg.db,
								 ont      = sub_ontology,
								 level    = 3,
								 readable = TRUE)
	if(results_only) return(get_GO_results(GO_object = ggo,sort_output = sort_results, sort_by = sort_results_by)) else return(ggo)
}


custom_enrichGO <- function(gene_list, sub_ontology="BP", results_only=FALSE, sort_results=FALSE, sort_results_by="Count"){
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
	if(results_only) return(get_GO_results(GO_object = ego,sort_output = sort_results, sort_by = sort_results_by)) else return(ego)
}


custom_enrichKEGG <- function(gene_list, results_only=FALSE, sort_results=FALSE, sort_results_by="Count"){
	# gene_list : vector of gene Entrez ID's
	kk <- enrichKEGG(gene         = gene, 
									 organism     = 'hsa',
									 pvalueCutoff = 0.05)
	if(results_only) return(get_GO_results(GO_object = kk,sort_output = sort_results, sort_by = sort_results_by)) else return(kk)
}


custom_gseKEGG <- function(gene_foldchange_list,results_only=FALSE, sort_results=FALSE, sort_results_by="setSize"){
	# gene_foldchange_list : vector of Fold Changes named by gene ID
	kk2 <- gseKEGG(geneList     = gene_foldchange_list, 
								 organism     = 'hsa',
								 nPerm        = 1000,
								 # minGSSize    = 120,
								 pvalueCutoff = 0.05,
								 verbose      = FALSE)
	if(results_only) return(get_GO_results(GO_object = kk2,sort_output = sort_results, sort_by = sort_results_by)) else return(kk2)
}


# extract the 'results' from the standard clusterProfiler output object
get_GO_results <- function(GO_object,sort_output=TRUE,sort_by="Count"){
	# save the results to a new dataframe
	ggo_results <- GO_object@result
	
	# sort by count
	if(sort_output) ggo_results <- ggo_results[with(ggo_results, order(get(sort_by),decreasing = TRUE)),] 
	return(ggo_results)
}



# custom barplot based on results dataframe
custom_GO_barplot <- function(GO_counts, GO_descriptions, plot_title="", plot_subtitle="", plot_margins=c(5, 4, 4, 2) + 0.1){
	# make a vector of values to plot
	GO_results_vec <- setNames(GO_counts,GO_descriptions)
	
	# set up colors for the bars in the plot
	plot_cols <- rainbow(n = length(GO_results_vec))
	
	# adjust the margins, these plots need lots of space
	# plot margins # c(bottom, left, top, right) # default is c(5, 4, 4, 2) + 0.1
	par(mar=plot_margins)
	barplot(GO_results_vec,horiz = TRUE,las=1,col = plot_cols,xlab = "Number of genes",cex.names=0.75,main = plot_title,sub = plot_subtitle)
}




# ~~~~~ GO CUSTOM WORKFLOW ~~~~~~~~~ #
out_dir <- "/ifs/home/kellys04/projects/SmithLab_ChIPSeq_2016-12-31/project_notes/gene_ontology"
diff_input_dir <- "/ifs/home/kellys04/projects/SmithLab_ChIPSeq_2016-12-31/project_notes/diffbind_files"

# element: file path, element name: Sample ID
diffbind_samples <- setNames(c(paste0(diff_input_dir,"/diff_bind.Sample4-ChIPSeq-vs-shControl-ChIPSeq.p100.csv"),
															 paste0(diff_input_dir,"/diff_bind.Sample4-ChIPSeq-vs-Sample5-ChIPSeq.p100.csv"),
															 paste0(diff_input_dir,"/diff_bind.Sample5-ChIPSeq-vs-shControl-ChIPSeq.p100.csv")),
														 c("Sample4-ChIPSeq-vs-shControl-ChIPSeq", 
														 	"Sample4-ChIPSeq-vs-Sample5-ChIPSeq",
														 	"Sample5-ChIPSeq-vs-shControl-ChIPSeq"))



# FULL WORKFLOW
# save all plots to a single PDF
# only process files 1, 3
pdf(file = paste0(out_dir,"/GO_plots6.pdf"),width = 10,height = 10,onefile = TRUE)
for(i in c(1,3)){
	tmp_sample <- diffbind_samples[i]
	message("Sample File:")
	message(tmp_sample)
	message("Sample Name:")
	message(names(tmp_sample))
	
	
	# read in the diffbind file
	message("Reading in DiffBind file")
	diffbind_df <- read.delim(file = tmp_sample,header = TRUE,sep = ',')
	# get the fold change values
	message("Getting Fold Change values")
	diffbind_foldchange <- get_diffbind_gene_foldchange(diffbind_df,sort_decreasing = TRUE)
	
	
	# run KEGG enrichment analysis
	message("Running KEGG Enrichment Analysis")
	diffbind_kegg_enrich_results <- custom_enrichKEGG(names(diffbind_foldchange),results_only = TRUE,sort_results = TRUE,sort_results_by = "Count")
	if(nrow(diffbind_kegg_enrich_results>=1)){
		message("Generating KEGG Enrichment Barplot")
		print(custom_GO_barplot(diffbind_kegg_enrich_results[["Count"]][1:10],diffbind_kegg_enrich_results[["Description"]][1:10], plot_title="KEGG Enrichment Plot",plot_subtitle=names(tmp_sample),plot_margins=c(5, 15, 4, 2) + 0.1))
	} 
	if(nrow(diffbind_kegg_enrich_results>=1)){
		message("Saving KEGG Enrichment table")
		write.table(x = diffbind_kegg_enrich_results, file = paste0(out_dir,"/",names(tmp_sample),"KEGG_enrichment.tsv"),quote = FALSE,sep = '\t',row.names = FALSE,col.names = TRUE)
	}  
	
	
	# run KEGG gene set enrichment analysis
	message("Running KEGG Gene Set Enrichment Analysis")
	diffbind_kegg_gse_results <- custom_gseKEGG(diffbind_foldchange,results_only = TRUE,sort_results = TRUE)
	if(nrow(diffbind_kegg_gse_results)>=1){
		message("Generating KEGG Gene Set Enrichment Barplot")
		print(custom_GO_barplot(diffbind_kegg_gse_results[["setSize"]][1:10],diffbind_kegg_gse_results[["Description"]][1:10], plot_title="KEGG Gene Set Enrichment Plot",plot_subtitle=names(tmp_sample),plot_margins=c(5, 15, 4, 2) + 0.1))
	}  
	if(nrow(diffbind_kegg_gse_results)>=1){
		message("Saving KEGG Gene Set Enrichment table")
		write.table(x = diffbind_kegg_gse_results, file = paste0(out_dir,"/",names(tmp_sample),"KEGG_gene_set_enrichment.tsv"),quote = FALSE,sep = '\t',row.names = FALSE,col.names = TRUE)
	}
	
	
	
	# run group GO functions for each subontology category
	ontologies <- setNames(c("MF","BP","CC"),c("Molecular Function", "Biological Process", "Cellular Component"))
	for(i in seq_along(ontologies)){
		sub_ont <- ontologies[i]
		# CC : cellular component
		# BP : biologacl process
		# MF : molecular funciton
		
		# run group GO functions
		message("Running group Gene Ontology analysis for: ")
		message(names(sub_ont))
		diffbind_group_gene_ontology_results <- custom_groupGO(names(diffbind_foldchange), sub_ont, results_only = TRUE,sort_results = TRUE,sort_results_by = "Count")
		if(nrow(diffbind_group_gene_ontology_results)>=1){
			message("Generating Group Gene Ontology barplot")
			print(custom_GO_barplot(diffbind_group_gene_ontology_results[["Count"]][1:10],diffbind_group_gene_ontology_results[["Description"]][1:10],plot_title = paste0("Group Gene Ontology plot: ",names(sub_ont)),plot_subtitle=names(tmp_sample), plot_margins = c(5, 20, 4, 2) + 0.1))
		}
		if(nrow(diffbind_group_gene_ontology_results)>=1){
			message("Saving Group Gene Ontology table")
			write.table(x = diffbind_group_gene_ontology_results, file = paste0(out_dir,"/",names(tmp_sample),"_",sub_ont,"_ggo.tsv"),quote = FALSE,sep = '\t',row.names = FALSE,col.names = TRUE)
		}
		
		
		# run enrichment GO functions
		message("Running Enrichment Gene Ontology analysis for: ")
		message(names(sub_ont))
		
		diffbind_enrich_gene_ontology_results <- custom_enrichGO(names(diffbind_foldchange), sub_ont, results_only = TRUE,sort_results = TRUE,sort_results_by = "Count")
		if(nrow(diffbind_enrich_gene_ontology_results)>=1){
			message("Generating Enrichment Gene Ontology barplot")
			print(custom_GO_barplot(diffbind_enrich_gene_ontology_results[["Count"]][1:10],diffbind_enrich_gene_ontology_results[["Description"]][1:10],plot_title = paste0("Enrichment Gene plot: ", names(sub_ont)),plot_subtitle=names(tmp_sample), plot_margins = c(5, 23, 4, 2) + 0.1))
		}
		if(nrow(diffbind_enrich_gene_ontology_results)>=1){
			message("Saving Enrichment Gene Ontology table")
			write.table(x = diffbind_enrich_gene_ontology_results, file = paste0(out_dir,"/",names(tmp_sample),"_",sub_ont,"_ego.tsv"),quote = FALSE,sep = '\t',row.names = FALSE,col.names = TRUE)
		}
		
		
		
	}
	
	
}
dev.off()

