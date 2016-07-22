#!/usr/bin/env Rscript

# quick script to make Volcano plots out of ChIP-Seq DiffBind data

# in the terminal you can use this command to get the columns:
# cat diff_bind.Sample4-ChIPSeq-vs-Control-ChIPSeq.p005.csv | sed 's/"//g' | awk '{FS=","; print $1":"$2"-"$3,$7,$8,$9,$10,$11}' > Volcano.xls 

# these are the columns we want:
# seqnames,start,end,width,strand,Conc,Conc_Sample4.ChIPSeq,Conc_Control.ChIPSeq,Fold,p.value,FDR,CBS1.Sample4.ChIPSeq,CBS2.Sample4.ChIPSeq,CBS3.Sample4.ChIPSeq,CBS7.Control.ChIPSeq,CBS8.Control.ChIPSeq,CBS9.Control.ChIPSeq,feature,external_gene_name,gene_biotype,start_position,end_position,insideFeature,distancetoFeature,shortestDistance,fromOverlappingOrNearest:-
# 	chr1:8158035-8158269 -0.31 3.82 -4.13 4.86e-06 1.04e-05
# chr1:10239463-10239824 4.31 1.9 2.42 3.91e-06 8.71e-06
# chr1:11346504-11346677 3.9 -0.37 4.27 2.26e-06 5.5e-06
# chr1:11346936-11347192 4.17 2.28 1.89 0.000155 0.000216
# chr1:11347260-11347505 3.92 1.18 2.74 5.58e-05 8.59e-05
# chr1:11415300-11415661 4.62 1.19 3.42 2.3e-09 1.39e-08
# chr1:11416116-11416402 3.93 -0.37 4.3 1.73e-06 4.36e-06
# chr1:11559582-11559821 4.07 0.7 3.37 2.68e-06 6.36e-06
# chr1:11762082-11762247 4.25 1.21 3.04 1.42e-06 3.68e-06



# place to save the plots
out_dir <- "/ifs/home/kellys04/projects/SmithLab_ChIPSeq-RNASeq_2016-12-31/project_notes/Combined_ChIP_RNA_Seq/analysis_output"


# place to get the DiffBind data from
pipeline_diffbind_results_dir <- "/ifs/home/kellys04/projects/SmithLab_ChIPSeq-RNASeq_2016-12-31/pipeline/diffbind/results/diffbind.by_chip.standard/peaks.by_sample.macs_narrow/align.by_sample.bowtie2/H33"

# dir(pipeline_diffbind_results_dir,pattern = "diff_bind*")
# > dir(pipeline_diffbind_results_dir,pattern = "diff_bind*")
# [1] "diff_bind.Sample4-ChIPSeq-vs-Control-ChIPSeq.p005.csv"     "diff_bind.Sample4-ChIPSeq-vs-Control-ChIPSeq.p020.csv"    
# [3] "diff_bind.Sample4-ChIPSeq-vs-Control-ChIPSeq.p100.csv"     "diff_bind.Sample4-ChIPSeq-vs-Sample5-ChIPSeq.p020.csv"
# [5] "diff_bind.Sample4-ChIPSeq-vs-Sample5-ChIPSeq.p100.csv" "diff_bind.Sample5-ChIPSeq-vs-Control-ChIPSeq.p005.csv"    
# [7] "diff_bind.Sample5-ChIPSeq-vs-Control-ChIPSeq.p020.csv"     "diff_bind.Sample5-ChIPSeq-vs-Control-ChIPSeq.p100.csv" 


# just work on one file for right now
diff_file <- "diff_bind.Sample4-ChIPSeq-vs-Sample5-ChIPSeq.p100.csv"


# read in the file
diff_df <- read.delim(file = paste0(pipeline_diffbind_results_dir,"/",diff_file),header = TRUE,sep = ",")


# check to colnames
colnames(diff_df)
# [1] "seqnames"                    "start"                       "end"                         "width"                      
# [5] "strand"                      "Conc"                        "Conc_Sample4.ChIPSeq"  "Conc_Control.ChIPSeq"     
# [9] "Fold"                        "p.value"                     "FDR"                         "CBS1.Sample4.ChIPSeq"
# [13] "CBS2.Sample4.ChIPSeq" "CBS3.Sample4.ChIPSeq" "CBS7.Control.ChIPSeq"     "CBS8.Control.ChIPSeq"    
# [17] "CBS9.Control.ChIPSeq"     "feature"                     "external_gene_name"          "gene_biotype"               
# [21] "start_position"              "end_position"                "insideFeature"               "distancetoFeature"          
# [25] "shortestDistance"            "fromOverlappingOrNearest"   


# make a subset of just the columns we want
plot_df <- diff_df[,c("seqnames","start","end","Conc_Sample4.ChIPSeq","Conc_Sample5.ChIPSeq","Fold","p.value","FDR")]

# NOTE:
# fold change = experimental / control
# only want pos. peaks; positive 2x Fold Change = log2(1) ; x values >=1 == 2x fold change

# subset again for just the values we want
plot_df_subset <- subset(plot_df, Fold >= 1)

# plot both to see how they look
plot(x = plot_df[["Fold"]],y = -log10(plot_df[["p.value"]]),pch=20, cex =0.5, col="blue", main="Volcano plot",xlab="log2(fold change)", xlim=c(-15,15),ylim=c(0,50))

plot(x = plot_df_subset[["Fold"]],y = -log10(plot_df_subset[["p.value"]]),pch=20, cex =0.5, col="blue", main="Volcano plot",xlab="log2(fold change)", xlim=c(-15,15),ylim=c(0,50))




