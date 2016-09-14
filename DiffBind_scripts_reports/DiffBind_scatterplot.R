#!/usr/bin/env Rscript

# create a scatter plot from DiffBind results, comparing differential binding between two samples
# x-axis: WT differential binding, y-axis: KO differential binding

DiffBind_scatter_plot_top_protein_coding_promoters <- function(diff_df, subset_FC = TRUE, subset_p = TRUE, subset_proteinCoding = TRUE, subset_promoter = TRUE, highlight_top_UP = TRUE, highlight_top_DOWN = TRUE, label_top_UP = TRUE, label_top_DOWN = TRUE, signif_p_value=0.05, fold_change_value=1, plot_colors=c("grey","red","orange","blue", "green"), top_gene_number=25){
    # get the colname order for later
    original_colnames <- colnames(diff_df)
    
    library("data.table")
    
    # make it a data.table type
    setDT(diff_df)
    # get min value per factor, also remove duplicates
    if(subset_promoter) diff_df <- as.data.frame(diff_df[, .SD[which.min(shortestDistance)], by=external_gene_name])
    # fix the colname order
    diff_df <- diff_df[,original_colnames]
    
    # get sample ID's
    sampleID_1 <- colnames(diff_df)[7]
    sampleID_2 <- colnames(diff_df)[8]
    print(sampleID_1)
    print(sampleID_2)
    
    # diff_df <- as.data.frame(subset(diff_df, p.value<signif_p_value & abs(Fold)>fold_change_value))
    
    # subset for significant genes
    if(subset_p) diff_df <- as.data.frame(subset(diff_df, p.value<signif_p_value))
    
    # subset for fold change           
    if(subset_FC) diff_df <-  as.data.frame(subset(diff_df, abs(Fold)>fold_change_value))
    
    # levels(diff_df[["gene_biotype"]])
    # "protein_coding"
    
    # subset for protein coding genes
    if(subset_proteinCoding) diff_df <- as.data.frame(subset(diff_df, gene_biotype=="protein_coding"))
    
    # subset top up/down protein coding genes
    diff_df_UpFold <- diff_df[with(diff_df, order(Fold, p.value)), c(sampleID_1,sampleID_2,"Fold","p.value","external_gene_name")][1:top_gene_number,]
    
    diff_df_DnFold <- diff_df[with(diff_df, order(-Fold, p.value)), c(sampleID_1,sampleID_2,"Fold","p.value","external_gene_name")][1:top_gene_number,]
    
    
    library("calibrate")
    # plot_lim <- max(na.rm = TRUE,abs(c(na.omit(diff_df[[sampleID_2]]), na.omit(diff_df[[sampleID_1]])))) + 2
    plot_lim <- 10
    print(plot_lim)
    print(plot_lim)
    # base plot, all points
    plot(x = diff_df[[sampleID_2]],
         y=diff_df[[sampleID_1]],
         xlab = sampleID_2,
         ylab = sampleID_1,
         col="grey",
         xlim = c(-1,plot_lim),
         ylim = c(-1,plot_lim),
         main = paste0(sampleID_1, " vs. ",sampleID_2))
    
    if(highlight_top_UP) {
        points(x = diff_df_UpFold[[sampleID_2]],
               y = diff_df_UpFold[[sampleID_1]],
               col='blue')
    }
    
    if(label_top_UP){
        textxy(diff_df_UpFold[[sampleID_2]],
               diff_df_UpFold[[sampleID_1]],
               labs=diff_df_UpFold[["external_gene_name"]])    
    }
    
    if(highlight_top_DOWN) {
        points(x = diff_df_DnFold[[sampleID_2]],
               y = diff_df_DnFold[[sampleID_1]],
               col='blue')
    }
    
    if(label_top_UP){
        textxy(diff_df_DnFold[[sampleID_2]],
               diff_df_DnFold[[sampleID_1]],
               labs=diff_df_DnFold[["external_gene_name"]])
    }
    
    abline(a = c(0, 1))
    abline(a = c(-1, 1), lty = 'dotted')
    abline(a = c(1, 1), lty = 'dotted')
    
    
    
    
    
}


out_dir <- "/ifs/home/kellys04/projects/SmithLab_ChIPSeq-RNASeq_2016-06-09/project_notes/diffbind_files"

diff_cl4_ctrl_file <- paste0(out_dir,"/diff_bind.Sample4-ChIPSeq-vs-Control-ChIPSeq.p100.csv")
diff_cl4_cl5_file <- paste0(out_dir,"/diff_bind.Sample4-ChIPSeq-vs-Sample5-ChIPSeq.p100.csv")
diff_cl5_ctrl_file <- paste0(out_dir,"/diff_bind.Sample5-ChIPSeq-vs-Control-ChIPSeq.p100.csv")

sample_files <- setNames(c(diff_cl4_ctrl_file, diff_cl4_cl5_file, diff_cl5_ctrl_file),
                         gsub(pattern = '.csv',replacement = '',x = gsub(pattern = 'diff_bind.',replacement = '',x = basename(c(diff_cl4_ctrl_file, diff_cl4_cl5_file, diff_cl5_ctrl_file)))))



# seq_along(names(sample_files)[c(1,3)])
for(i in c(1,3)){
    diff_file <- sample_files[i]
    print(diff_file)
    
    diff_df <- read.delim(file = diff_file, header = TRUE,sep = ',')
    
    outfile_path_all <- paste0(out_dir,"/scatterPlot_",names(sample_files[i]),"_all-ProteinPromoters.pdf")
    outfile_path_filtered <- paste0(out_dir,"/scatterPlot_",names(sample_files[i]),"_filtered-FC_P.pdf")
    
    pdf(file = outfile_path_filtered,height = 8,width = 8)
    print(DiffBind_scatter_plot_top_protein_coding_promoters(diff_df, label_top_UP = FALSE, label_top_DOWN = FALSE))
    dev.off()
    
    pdf(file = outfile_path_all,height = 8,width = 8)
    print(DiffBind_scatter_plot_top_protein_coding_promoters(diff_df,subset_FC = FALSE,subset_p = FALSE, label_top_UP = FALSE, label_top_DOWN = FALSE))
    dev.off()
    
}

