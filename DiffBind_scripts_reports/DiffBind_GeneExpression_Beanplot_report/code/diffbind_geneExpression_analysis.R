#!/usr/bin/env Rscript

## USAGE: diffbind_geneExpression_analysis.R /path/to/outdir /path/to/microarray_gene_expression.tsv /path/to/diffbind.csv <HistMark> <analysis_params_branch>

## OUTPUT: <outdir>/<SampleID>/<HistMark>/<analysis_params_branch>/[paired_boxplot.pdf genes_up.tsv genes_down.tsv]
## e.g. <outdir>/ZNK/INPUT/peaks.by_sample.macs_broad-cutoff0025

## DESCRIPTION:
# This pipeline script will compare microarray gene expression data with
# DiffBind differential peak binding data
# for each sample, the script will find the set of genes that are up or down regulated in the microarray data
# and create subsets of genes with differentially bound peaks from the DiffBind data
# and plot the DiffBind values for both subsets (up and down regulated genes) in a paired boxplot
# the subset diffbind values will be output as TSV, and the plot

# ~~~~~~ file dir location notes ~~~~~~~~ #
# on server:
# pipeline_diffbind_dir <- "/ifs/data/sequence/results/smithlab/2016-01-04/ChIP-Seq/diffbind/results/diffbind.by_chip.status/peaks.by_sample.macs_broad/align.by_sample.bowtie2"
# main_outdir <- "/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2016-03-10/project_notes/diffferential_gene_expression/diff_expr_report/diff_expr_TSS_EdgeR"
# diffbind_file <- "/ifs/data/sequence/results/smithlab/2016-01-04/ChIP-Seq/diffbind/results/diffbind.by_chip.status/peaks.by_sample.macs_broad/align.by_sample.bowtie2/H3K27AC/diff_bind.D-vs-R.p100.csv"
# # gene_expr_file <- "/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2016-03-10/project_notes/diffferential_gene_expression/GeneExpression18Patients_includes_SRR.txt"
# gene_expr_file <- "/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2016-03-10/project_notes/diffferential_gene_expression/gene_expression.tsv"

# local:
# pipeline_diffbind_dir <- "/ifs/data/sequence/results/smithlab/2016-01-04/ChIP-Seq/diffbind/results/diffbind.by_chip.status/peaks.by_sample.macs_broad/align.by_sample.bowtie2"
# # main_outdir <- "/Users/kellys04/projects/tmp_Smith/diffferential_gene_expression/diff_expr_TSS_EdgeR"
# setwd(main_outdir)
# diffbind_file <- "/Users/kellys04/projects/tmp_Smith/diffferential_gene_expression/diff_bind.D-vs-R.p100.csv"
# gene_expr_file <- "/Users/kellys04/projects/tmp_Smith/diffferential_gene_expression/gene_expression.tsv"
# histone_mark <- "H3K27AC"
# params_branch <- "peaks.by_sample.macs_broad"

# set up the analysis_output lke this:
# ZNK/INPUT/peaks.by_sample.macs_broad-cutoff0025
# ~~~~~~~~~~~~~~~ #



 
# install.packages(c("devtools"))
# source("http://www.bioconductor.org/biocLite.R")
# biocLite(c("Biobase","preprocessCore"))
library("data.table")
library("beanplot")
library("ggplot2")
library("reshape2")
library("preprocessCore")




# ~~~~~ CUSTOM FUNCTIONS ~~~~~~~~ #
get_replacements <- function(source_data, uniq=TRUE, pattern = "", remove_entry = NULL, ...){
    # return vector of first matches to regex
    matches <-  gsub(x = source_data, 
                     pattern = pattern,
                     replacement = "\\1")
    
    # return only unique entries
    if(uniq) matches <- unique(matches)
    
    # remove an entry from the output
    if ( ! is.null(remove_entry)){
        if(any(grepl(pattern = remove_entry,x = matches))){
            matches <- matches[grep(pattern = remove_entry, x = matches, value = FALSE, invert = TRUE)]  
        } 
    }
    return(matches)
}


multi_grep <- function(source_data, patterns){
    # find multiple patterns in a char vector
    # make a regex to search with
    if (length(patterns) > 1){
        search_pattern <- paste(patterns,
                                collapse = "|") 
    } else {
        search_pattern <- patterns
    }
    
    matches <- grep(x = source_data,
                    pattern = search_pattern,
                    value = TRUE)
    return(matches)
}

sampleID_intersects <- function(IDs1, pattern1 = "", IDs2, pattern2 = "", removeID2 = NULL ){
    # get the intersect between sets of ID's after gsub pattern replacements
    ID1_match <- get_replacements(source_data = IDs1, pattern = pattern1)
    
    ID2_intersect_ID1 <-  multi_grep(patterns = ID1_match, source_data = IDs2)
    ID2_match <- get_replacements(source_data = ID2_intersect_ID1, pattern = pattern2, remove_entry = removeID2)
    return(ID2_match)
}

quantile_normalize_df <- function(df){
    # quantile normalization of a dataframe
    df <- setnames(x = as.data.frame(normalize.quantiles(as.matrix(df))), colnames(df))
    return(df)
}

diff_subset <- function(df, patientID_matches){
    # subset for DiffPeaks >= 1 per sample
    for(i in seq_along(patientID_matches)){
        tmp_ID <- patientID_matches[i]
        
        # print(tmp_ID %in% colnames(diffbind_df_intersect_min))
        
        # get the diffBind columns with the patient's ID
        tmp_colnames <- grep(pattern = tmp_ID,x = colnames(df),value = TRUE)
        
        # make sure there are 2 colnames, otherwise skip to the next iteration of the loop!
        if (length(tmp_colnames) !=2) next
        
        # get the R colname
        tmp_colnames_R <- grep(pattern = paste0("^",tmp_ID,"[[:punct:]]R.*$"),x = tmp_colnames,perl = TRUE,value = TRUE)
        
        # get the D colname
        tmp_colnames_D <- grep(pattern = paste0("^",tmp_ID,"[[:punct:]]D.*$"),x = tmp_colnames,perl = TRUE,value = TRUE)
        
        # diffbind_df <- subset(diffbind_df, tmp_colnames_R >= 1)
        # diffbind_df <- subset(diffbind_df, tmp_colnames_D >= 1)
        df <- df[which(df[[tmp_colnames_R]] >= 1 | df[[tmp_colnames_D]] >= 1),]
        # diffbind_df <- diffbind_df[which(diffbind_df[[tmp_colnames_D]] >= 1),]
    }
    
    return(df)
}

quantile_cutoff <- function(vec, cutoff = 0.25){
    # replace lowest values with the x'th quantile value
    vec_quantile <- quantile(vec, probs = cutoff)
    vec <- sapply(vec, function(x) max(x, vec_quantile))
    return(vec)
}


calc_fold_change <- function(sample1, sample2, calc_log = TRUE, log_base = 2){
    # sample1, sample2: numeric vectors
    # fc: fold change output vector
    
    if(calc_log){
        fc <- log((sample2/sample1), base = log_base)  
    } else {
        fc <- (sample2/sample1)
    }
    
    return(fc)
}


log_ratio_df <- function(diffbind_df, diff_gene_colname = "external_gene_name", patientID_matches, calc_cutoff = FALSE, melt_df = FALSE, add_status = FALSE, ...){
    # create a df to hold the DiffBind log ratio entries
    # create an empty data frame to fill with log ratio values
    diff_log_ratio <- setNames(data.frame(matrix(nrow = nrow(diffbind_df), ncol = length(patientID_matches))), patientID_matches)
    # add the genes from DiffBind
    diff_log_ratio[["gene"]] <- diffbind_df[[diff_gene_colname]]
    
    # enter the log ratio values into the df per patient
    # # iterate over the patients; get only diff min cols per patient
    for(i in seq_along(patientID_matches)){
        tmp_ID <- patientID_matches[i]
        
        # get the diffBind columns with the patient's ID
        tmp_colnames <- grep(pattern = tmp_ID,x = colnames(diffbind_df),value = TRUE)
        
        # make sure there are 2 colnames, otherwise skip to the next iteration of the loop!
        if (length(tmp_colnames) !=2) next
        
        # get the R colname
        tmp_colnames_R <- grep(pattern = paste0("^",tmp_ID,"[[:punct:]]R.*$"),x = tmp_colnames,perl = TRUE,value = TRUE)
        
        # get the D colname
        tmp_colnames_D <- grep(pattern = paste0("^",tmp_ID,"[[:punct:]]D.*$"),x = tmp_colnames,perl = TRUE,value = TRUE)
        
        if(calc_cutoff){
            # apply quantil cutoffs to the R and D values
            tmp_R <- quantile_cutoff(diffbind_df[[tmp_colnames_R]])
            tmp_D <- quantile_cutoff(diffbind_df[[tmp_colnames_D]])
        } else {
            tmp_R <- diffbind_df[[tmp_colnames_R]]
            tmp_D <- diffbind_df[[tmp_colnames_D]]
        }
        
        # get the fold change
        diff_log_ratio[,tmp_ID] <- calc_fold_change(sample1 = tmp_D, 
                                                    sample2 = tmp_R, 
                                                    calc_log = TRUE, log_base = 2)
        
    }
    
    
    # melt it into long format
    if (melt_df) diff_log_ratio <- reshape2::melt(diff_log_ratio,id.vars="gene",variable.name="sample",value.name="diff_peak_log_ratio")
    
    # add a column to hold the 'status' e.g. Up regulated vs. Down regulated in the gene_expr table
    if (add_status) diff_log_ratio[["gene_expression_status"]] <- NA
    
    return(diff_log_ratio)
}

mark_test_type <- function(hist_mark, mark_key){
    # get the alternative value for the ttest/utest from the mark key
    test_alt <- as.character(subset(mark_key, marks == hist_mark)[['test']])
    return(test_alt)
}

plot_signif_pvalue <- function(df, y_colname, x_colname, x_group, p_value_cutoff=0.05, test_method = "ttest", ...){
    # plot points on a barplot / beanplot marking significant categories
    y_max <- max(df[[y_colname]],na.rm = TRUE)
    
    for(i in seq_along(unique(as.character(df[[x_group]])))){
        tmp_ID <- unique(as.character(df[[x_group]]))[i]
        tmp_df <- subset(df, get(x_group) == tmp_ID)
        
        cat(test_method)
        test_alt <- mark_test_type(...)
        cat(test_alt)
        if(test_method == 'ttest') tmp_pvalue <- t.test(get(y_colname)~get(x_colname),
                                                        data = tmp_df, 
                                                        alternative = test_alt,
                                                        na.action = na.omit,
                                                        paired = FALSE)["p.value"]
        
        # Wilcoxon rank-sum test applied to independent samples
        if(test_method == 'utest'){
            tmp_pvalue <- wilcox.test(get(y_colname)~get(x_colname),
                                      data = tmp_df,
                                      alternative = test_alt,
                                      na.action = na.omit, 
                                      paired = FALSE)["p.value"]
            
            print(wilcox.test(get(y_colname)~get(x_colname),
                              data = tmp_df,
                              paired = FALSE))
            
        } 
        if(tmp_pvalue < p_value_cutoff){
            print(points(x = i, y = y_max - 1,pch ='*', col ='red',cex=2))
            print(text(x = i, y = y_max,labels = paste0('p = ',format(tmp_pvalue, digits = 2)),cex=0.7))
        }
    }
}



diff_beanplot <- function(df, x1_colname, x2_colname, y_colname, main_text_1 = "", y_lab = "",
                          main_text_2 = "", save_output = FALSE, file = "./plot.pdf", strip_chart = FALSE, ... ){
    # make the DiffBind split beanplots
    
    # only make a plot of there are at least 50 genes
    if(nrow(unique(subset(df, get(x1_colname) == "UP")["gene"]))>=50 &&
       nrow(unique(subset(df, get(x1_colname) == "DOWN")["gene"]))>=50){
        
        if(save_output) pdf(file = file,height = 8,width = 12)
        
        beanplot(get(y_colname)~get(x1_colname)*get(x2_colname),
                 data=df,
                 what=c(0,1,1,0), # what=c(0,1,1,1),maxstripline = 0, # ll = 0.04,varwidth = TRUE,
                 border = NA,
                 bw="nrd0", # had to change this, was giving errors with default and nrd 
                 overallline = 'mean', #median',
                 col=list('grey','purple'),
                 ylab = y_lab, 
                 main=paste0("DiffBind results for ",main_text_1," ",main_text_2), cex.main=0.9
                 ,side = "both"
        )
        
        legend('bottomright', fill=c('grey','purple'), legend= c('Down Genes', 'Up Genes'))
        
        abline(h=0,col="darkgrey",lty="dotted")
        
        # add p value and *'s
        plot_signif_pvalue(df = df,y_colname = y_colname, x_colname = x1_colname, x_group = x2_colname, ...) # , test_method = 'utest'
        
        if (strip_chart){
            stripchart(y_colname~x_colname*sample,data = df,
                       vertical=TRUE,
                       add=TRUE,
                       cex=0.1,
                       method="jitter",
                       jitter=0.05)   
        }
        
        if(save_output) dev.off()
    }
}

gene_expression_raw_pipeline <- function(gene_expr_df, patientID_matches, gene_names_overlap, quant_norm = FALSE){
    # keep only LR columns; log ratio
    gene_expr_df <- gene_expr_df[multi_grep(source_data = colnames(gene_expr_df), patterns = "^.*LR$")]
    # quantil normalize the gene expression dataset
    if (quant_norm) gene_expr_df <- quantile_normalize_df(gene_expr_df)
    # fix colnames
    colnames(gene_expr_df) <- get_replacements(source_data = colnames(gene_expr_df), uniq=TRUE, pattern = "^([[:alpha:]]*)([[:punct:]]*.*)$")
    # keep matched columns
    gene_expr_df <- gene_expr_df[patientID_matches]
    # add a column with the gene names for melting
    gene_expr_df[["gene"]] <- rownames(gene_expr_df)
    
    # subset for intersected genes only
    gene_expr_df <- droplevels(gene_expr_df[rownames(gene_expr_df)  %in% gene_names_overlap,])
    
    # melt it into long format
    gene_expr_df <- reshape2::melt(gene_expr_df,id.vars="gene",variable.name="sample",value.name="gene_expression_log_ratio")
    
    return(gene_expr_df)
}

diffbind_raw_pipeline <- function(diffbind_df, patientID_matches, gene_names_overlap, quant_norm = FALSE){
    # quantile normalize entires DiffBind dataset
    if (quant_norm) {
        norm_cols <- multi_grep(source_data = colnames(diffbind_df), patterns = patientID_matches)
        diffbind_df[norm_cols] <- quantile_normalize_df(diffbind_df[norm_cols])
    }
    
    
    # subset for closest genes
    diffbind_df <- subset(diffbind_df, abs(distancetoFeature)<=3000)
    # subset for intersected genes only
    diffbind_df <- droplevels(diffbind_df[diffbind_df[["external_gene_name"]] %in% gene_names_overlap,])
    # remove extraneous columns
    diffbind_cols_to_remove <-c("seqnames","start","end","width","strand","Conc","Conc_D","Conc_R","Fold",
                                "p.value","FDR","feature","gene_biotype","start_position","end_position",
                                "insideFeature","shortestDistance","fromOverlappingOrNearest", "distancetoFeature")
    diffbind_df <- diffbind_df[,! colnames(diffbind_df) %in% diffbind_cols_to_remove]
    return(diffbind_df)
}



combined_full_pipeline <- function(diffbind_df, gene_expr_df, patientID_matches, gene_names_overlap, 
                                   diff_quant_norm = FALSE, gene_quant_norm = FALSE, diff_lr_cutoff = FALSE, 
                                   make_plot = TRUE, data_file = "./data.Rdata", plot_file = "./plot.pdf", histone_mark, 
                                   params_branch, mark_key, ...){
    # ~~~~~ MANIPULATE RAW DATA ~~~~~~~~ #
    # DIFFBIND DATA
    diffbind_df <- diffbind_raw_pipeline(diffbind_df = diffbind_df, 
                                         patientID_matches = patientID_matches, 
                                         gene_names_overlap = gene_names_overlap, 
                                         quant_norm = diff_quant_norm)
    # GENE EXPRESSION DATA
    gene_expr_df <- gene_expression_raw_pipeline(gene_expr_df = gene_expr_df, 
                                                 patientID_matches = patientID_matches,
                                                 gene_names_overlap = gene_names_overlap, 
                                                 quant_norm = gene_quant_norm)
    
    # ~~~~~ CALCULATE DIFFBIND LOG RATIO  ~~~~~~~~ #
    diff_log_ratio <- log_ratio_df(diffbind_df = diffbind_df, 
                                   diff_gene_colname = "external_gene_name", 
                                   patientID_matches = patientID_matches, 
                                   calc_cutoff = diff_lr_cutoff,
                                   melt_df = TRUE, add_status = TRUE)
    
    # ~~~~~ MERGE THE TABLES ~~~~~~~~ #
    # merge together the two long format diff_log_ratio's and gene_expr log ratio's
    diff_gene_expr_merged <- base::merge(gene_expr_df,diff_log_ratio,by=c("gene","sample")) # ,all=TRUE
    # set the gene expression status to UP or DOWN based on gene expression value
    # 1.5x up/down expression
    diff_gene_expr_merged[["gene_expression_status"]] <- ifelse(diff_gene_expr_merged[["gene_expression_log_ratio"]]>=0.58, "UP",
                                                                ifelse(diff_gene_expr_merged[["gene_expression_log_ratio"]]<=-0.58,"DOWN",no = NA))
    # ~~~~~~~~~~~~~~~ #
    
    # ~~~~~~ SUBSET DIFF FOLD CHANGE ~~~~~~~~~ #
    # subset for log ratio >1.5x change; log2 0.58
    # diff_gene_expr_merged <- subset(diff_gene_expr_merged, abs(diff_peak_log_ratio) >= 0.58)
    # diff_gene_expr_merged <- droplevels(diff_gene_expr_merged)
    # ~~~~~~~~~~~~~~~ #
    
    
    if (make_plot){
        
        # ~~~~~ PLOT DIFFBIND VALUES ~~~~~~~~ #
        save.image(file=data_file, compress = TRUE) # 
        
        
        diff_beanplot(df = diff_gene_expr_merged,
                      x1_colname =  'gene_expression_status', 
                      x2_colname = 'sample',
                      y_colname = 'diff_peak_log_ratio', 
                      main_text_1 = histone_mark, 
                      main_text_2 = params_branch,
                      y_lab = "DiffBind ratio = log2 ( Diff peak R / Diff peak D )",
                      test_method = 'utest',
                      hist_mark = histone_mark, 
                      mark_key = mark_key,
                      save_output = TRUE, file = plot_file) # file = 
        # ~~~~~~~~~~~~~~~ #
    }
}

# ~~~~~~~~~~~~~~~ #







# ~~~~~~ GET SCRIPT ARGS ~~~~~~~~ #
## USAGE: diffbind_geneExpression_analysis.R /path/to/outdir /path/to/microarray_gene_expression.tsv /path/to/diffbind.csv <HistMark> <analysis_params_branch>
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

# histone marks included in the experiment
mark_key <- data.frame(marks = c('H3K9ME3', 'H3K27ME3', 'H3K9AC', 'H3K27AC', 'H3K4ME3'), 
                       type = c('repressive', 'repressive','activating','activating','activating'),
                       test = c('greater', 'greater', 'less', 'less', 'less'))
# ~~~~~~~~~~~~~~~ #


# ~~~~~~ READ IN FILE ~~~~~~~~ #
# diffbind data
diffbind_df <- read.csv(diffbind_file)
# gene expression data
gene_expr_df <- read.table(gene_expr_file,sep = '\t',header = TRUE,row.names = 1,quote = "")

# save current state
save.image(file=paste0(main_outdir,"/raw_data.Rdata"),compress = TRUE)
# ~~~~~~~~~~~~~~~ #


# ~~~~~ GET DATA INTERSECTS ~~~~~~~~ #
patientID_matches <- sampleID_intersects(IDs1 = colnames(gene_expr_df), 
                                         pattern1 = "^([[:alpha:]]*)[[:punct:]]*.*$", 
                                         IDs2 = colnames(diffbind_df), 
                                         pattern2 = "^([[:alpha:]]*)[[:punct:]]*.*$", 
                                         removeID2 = "GHW")

# get the genes in common between the diffbind and gene expression data
gene_names_overlap <- intersect(rownames(gene_expr_df),unique(diffbind_df[["external_gene_name"]]))
# ~~~~~~~~~~~~~~~ #



# without normalization
combined_full_pipeline(diffbind_df = diffbind_df, 
                       gene_expr_df = gene_expr_df, 
                       patientID_matches = patientID_matches, 
                       gene_names_overlap = gene_names_overlap,
                       diff_quant_norm = FALSE, 
                       gene_quant_norm = FALSE, 
                       diff_lr_cutoff = FALSE,
                       make_plot = TRUE, 
                       data_file = paste0(main_outdir,"/data.Rdata"), 
                       plot_file = paste0(main_outdir,"/beanplot_no_normalization.pdf"), 
                       histone_mark = histone_mark, 
                       params_branch = params_branch, 
                       mark_key = mark_key)

# DiffBind matrix quantile normalization
combined_full_pipeline(diffbind_df = diffbind_df, 
                       gene_expr_df = gene_expr_df, 
                       patientID_matches = patientID_matches, 
                       gene_names_overlap = gene_names_overlap,
                       diff_quant_norm = TRUE, 
                       gene_quant_norm = FALSE, 
                       diff_lr_cutoff = FALSE,
                       make_plot = TRUE, 
                       data_file = paste0(main_outdir,"/data.Rdata"), 
                       plot_file = paste0(main_outdir,"/beanplot_DiffQnorm.pdf"), 
                       histone_mark = histone_mark, 
                       params_branch = params_branch, 
                       mark_key = mark_key)

# DiffBind matrix quantile normalization + per sample quantile cutoffs
combined_full_pipeline(diffbind_df = diffbind_df, 
                       gene_expr_df = gene_expr_df, 
                       patientID_matches = patientID_matches, 
                       gene_names_overlap = gene_names_overlap,
                       diff_quant_norm = TRUE, 
                       gene_quant_norm = FALSE, 
                       diff_lr_cutoff = TRUE,
                       make_plot = TRUE, 
                       data_file = paste0(main_outdir,"/data.Rdata"), 
                       plot_file = paste0(main_outdir,"/beanplot_DiffQnorm_Qcutoff.pdf"), 
                       histone_mark = histone_mark, 
                       params_branch = params_branch, 
                       mark_key = mark_key)

# per sample quantile cutoffs
combined_full_pipeline(diffbind_df = diffbind_df, 
                       gene_expr_df = gene_expr_df, 
                       patientID_matches = patientID_matches, 
                       gene_names_overlap = gene_names_overlap,
                       diff_quant_norm = FALSE, 
                       gene_quant_norm = FALSE, 
                       diff_lr_cutoff = TRUE,
                       make_plot = TRUE, 
                       data_file = paste0(main_outdir,"/data.Rdata"), 
                       plot_file = paste0(main_outdir,"/beanplot_Qcutoff.pdf"), 
                       histone_mark = histone_mark, 
                       params_branch = params_branch, 
                       mark_key = mark_key)


# DiffBind matrix quantile normalization + per sample quantile cutoffs + gene expression matrix quantile normalization
combined_full_pipeline(diffbind_df = diffbind_df, 
                       gene_expr_df = gene_expr_df, 
                       patientID_matches = patientID_matches, 
                       gene_names_overlap = gene_names_overlap,
                       diff_quant_norm = TRUE, 
                       gene_quant_norm = TRUE, 
                       diff_lr_cutoff = TRUE,
                       make_plot = TRUE, 
                       data_file = paste0(main_outdir,"/data.Rdata"), 
                       plot_file = paste0(main_outdir,"/beanplot_DiffQnorm_GExprQnorm_Qcutoff.pdf"), 
                       histone_mark = histone_mark, 
                       params_branch = params_branch, 
                       mark_key = mark_key)



