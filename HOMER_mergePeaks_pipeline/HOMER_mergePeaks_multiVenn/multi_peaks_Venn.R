#!/usr/bin/env Rscript

## USAGE: multi_peaks_Venn.R <sampleID> /path/to/venn.txt
## This script will process venn output from HOMER mergePeaks and plot a venn diagram
## this is currently only for R version 3.2.0; module load r/3.2.0

# ~~~~~ LOAD PACKAGES ~~~~~~~ #
library('VennDiagram')
library('gridExtra')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 


# ~~~~~ GET SCRIPT ARGS ~~~~~~~ #
args <- commandArgs(TRUE); cat("Script args are:\n"); args
SampleID<-args[1]
venn_table_file<-args[2]
plot_outdir<-dirname(venn_table_file)
plot_filepath<-paste0(plot_outdir,"/",SampleID,"_venn.pdf") 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 


# ~~~~~ PARSE THE VENN TABLE ~~~~~~~ #
# read in the venn text
venn_table_df<-read.table(venn_table_file,header = TRUE,sep = "\t",stringsAsFactors = FALSE,check.names = FALSE)
# venn_table_df

# get the venn categories
venn_categories<-colnames(venn_table_df)[!colnames(venn_table_df) %in% c("Total","Name")] 
cat("Venn categories are:\n"); venn_categories

# venn_categories
num_categories<-length(venn_categories)
cat("Num categories are:\n"); num_categories

# make a summary table
venn_summary<-venn_table_df[!colnames(venn_table_df) %in% venn_categories]
cat("Venn summary table is categories are:\n"); venn_summary
# venn_summary

# write summary table
write.table(venn_summary,file = "venn_summary.tsv",quote = FALSE,row.names = FALSE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 


# ~~~~~ SET UP THE PLOT ~~~~~~~ #
# get the areas for the venn; add up all the overlaps that contain the given category 

if (num_categories == 2) {
  # PAIRWISE VENN
  cat("CREATING PAIR-WISE VENN DIAGRAM\n")
  # area1
  area_n1<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")"),
                                 x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # area2
  area_n2<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[2],")"),
                                 x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n12
  area_n12<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[2],")"),
                                  x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  venn <-draw.pairwise.venn(area1=area_n1,
                          area2=area_n2,
                          cross.area=area_n12,
                          category=gsub(pattern = ".bed",replacement = "",x = venn_categories),
                          fill=c('red','blue'),
                          alpha=0.3,
                          # cat.dist = 0.1,
                          cex=2,
                          cat.cex = 2,
                          margin = 0.1,
                          ind = FALSE)
  
  pdf(plot_filepath,width = 8,height = 8)
  # grid.arrange(gTree(children=venn), top=paste0(SampleID," Peak Overlap")) #, bottom="subtitle")
  grid.arrange(gTree(children=venn), top=textGrob(paste0(SampleID," Peak Overlap"), gp=gpar(cex=2), vjust=3))
  # top=textGrob("Total Data and Image", gp=gpar(cex=3), just="top") #  vjust=0.7
  dev.off()
  
  
} else if ( num_categories == 3 ) {
  # 3-WAY VENN
  cat("CREATING THREE-WAY VENN DIAGRAM\n")
  
  # area1
  area_n1<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")"),
                                 x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # area2
  area_n2<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[2],")"),
                                 x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # area3
  area_n3<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[3],")"),
                                 x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n12
  area_n12<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[2],")"),
                                  x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n13
  area_n13<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[3],")"),
                                  x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n23
  area_n23<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[2],")","(?=.*",venn_categories[3],")"),
                                  x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  
  # n123
  area_n123<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[2],")","(?=.*",venn_categories[3],")"),
                                   x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  
  venn <-draw.triple.venn(area1=area_n1,
                          area2=area_n2,
                          area3=area_n3,
                          n12=area_n12,
                          n13=area_n13,
                          n23=area_n23,
                          n123=area_n123,
                          category=gsub(pattern = ".bed",replacement = "",x = venn_categories),
                          fill=c('red','blue','green'),
                          alpha=0.3,
                          cex=2,
                          ind = FALSE)
  
  pdf(plot_filepath,width = 9,height = 9)
  grid.arrange(gTree(children=venn), top=paste0(SampleID," Peak Overlap")) #, bottom="subtitle")
  dev.off()
  
} else if ( num_categories == 4 ) {
  # 4-WAY VENN
  cat("CREATING FOUR-WAY VENN DIAGRAM\n")
  # area1
  area_n1<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")"),
                                 x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # area2
  area_n2<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[2],")"),
                                 x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # area3
  area_n3<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[3],")"),
                                 x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # area4
  area_n4<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[4],")"),
                                 x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  
  # n12
  area_n12<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[2],")"),
                                  x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n13
  area_n13<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[3],")"),
                                  x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n14
  area_n14<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[4],")"),
                                  x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  
  # n23
  area_n23<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[2],")","(?=.*",venn_categories[3],")"),
                                  x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n24
  area_n24<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[2],")","(?=.*",venn_categories[4],")"),
                                  x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  
  # n34
  area_n34<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[3],")","(?=.*",venn_categories[4],")"),
                                  x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  
  
  # n123
  area_n123<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[2],")","(?=.*",venn_categories[3],")"),
                                   x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n124
  area_n124<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[2],")","(?=.*",venn_categories[4],")"),
                                   x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  
  # n134
  area_n134<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[3],")","(?=.*",venn_categories[4],")"),
                                   x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n234
  area_n234<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[2],")","(?=.*",venn_categories[3],")","(?=.*",venn_categories[4],")"),
                                   x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n1234
  area_n1234<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[2],")","(?=.*",venn_categories[3],")",
                                                     "(?=.*",venn_categories[4],")"),
                                    x = venn_summary$Name,perl = TRUE),][["Total"]])
  

  venn <- draw.quad.venn(
    area1 = area_n1,
    area2 = area_n2,
    area3 = area_n3,
    area4 = area_n4,
    n12 = area_n12,
    n13 = area_n13,
    n14 = area_n14,
    n23 = area_n23,
    n24 = area_n24,
    n34 = area_n34,
    n123 = area_n123,
    n124 = area_n124,
    n134 = area_n134,
    n234 = area_n234,
    n1234 = area_n1234,
    category = gsub(pattern = ".bed",replacement = "",x = venn_categories),
    fill = c("orange", "red", "green", "blue"),
    cat.dist = 0.25,
    cat.cex = 1.2,
    alpha=0.3,
    margin = 0.1,
    ind = FALSE)

  pdf(plot_filepath,width = 9,height = 9)
  grid.arrange(gTree(children=venn), top=paste0(SampleID," Peak Overlap")) #, bottom="subtitle")
  dev.off()
  
} else if ( num_categories == 5 ) {
  # 5-WAY VENN
  cat("CREATING FIVE-WAY VENN DIAGRAM\n")
  
  # area1
  area_n1<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")"),
                                 x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # area2
  area_n2<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[2],")"),
                                 x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # area3
  area_n3<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[3],")"),
                                 x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # area4
  area_n4<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[4],")"),
                                 x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # area5
  area_n5<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[5],")"),
                                 x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  
  
  # n12
  area_n12<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[2],")"),
                                  x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n13
  area_n13<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[3],")"),
                                  x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n14
  area_n14<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[4],")"),
                                  x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n15
  area_n15<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[5],")"),
                                  x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n23
  area_n23<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[2],")","(?=.*",venn_categories[3],")"),
                                  x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n24
  area_n24<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[2],")","(?=.*",venn_categories[4],")"),
                                  x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n25
  area_n25<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[2],")","(?=.*",venn_categories[5],")"),
                                  x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n34
  area_n34<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[3],")","(?=.*",venn_categories[4],")"),
                                  x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n35
  area_n35<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[3],")","(?=.*",venn_categories[5],")"),
                                  x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n45
  area_n45<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[4],")","(?=.*",venn_categories[5],")"),
                                  x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  
  
  
  # n123
  area_n123<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[2],")","(?=.*",venn_categories[3],")"),
                                   x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n124
  area_n124<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[2],")","(?=.*",venn_categories[4],")"),
                                   x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n125
  area_n125<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[2],")","(?=.*",venn_categories[5],")"),
                                   x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n134
  area_n134<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[3],")","(?=.*",venn_categories[4],")"),
                                   x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n135
  area_n135<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[3],")","(?=.*",venn_categories[5],")"),
                                   x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n145
  area_n145<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[4],")","(?=.*",venn_categories[5],")"),
                                   x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  
  # n234
  area_n234<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[2],")","(?=.*",venn_categories[3],")","(?=.*",venn_categories[4],")"),
                                   x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n235
  area_n235<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[2],")","(?=.*",venn_categories[3],")","(?=.*",venn_categories[5],")"),
                                   x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n245
  area_n245<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[2],")","(?=.*",venn_categories[4],")","(?=.*",venn_categories[5],")"),
                                   x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n345
  area_n345<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[3],")","(?=.*",venn_categories[4],")","(?=.*",venn_categories[5],")"),
                                   x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  
  
  # n1234
  area_n1234<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[2],")","(?=.*",venn_categories[3],")",
                                                     "(?=.*",venn_categories[4],")"),
                                   x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n1235
  area_n1235<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[2],")","(?=.*",venn_categories[3],")",
                                                     "(?=.*",venn_categories[5],")"),
                                    x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n1245
  area_n1245<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[2],")","(?=.*",venn_categories[4],")",
                                                     "(?=.*",venn_categories[5],")"),
                                    x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  # n1345
  area_n1345<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[3],")","(?=.*",venn_categories[4],")",
                                                     "(?=.*",venn_categories[5],")"),
                                    x = venn_summary$Name,perl = TRUE),][["Total"]])
  # n2345
  area_n2345<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[2],")","(?=.*",venn_categories[3],")","(?=.*",venn_categories[4],")",
                                                     "(?=.*",venn_categories[5],")"),
                                    x = venn_summary$Name,perl = TRUE),][["Total"]])
  
  
  # n12345
  area_n12345<-sum(venn_summary[grep(pattern = paste0("(?=.*",venn_categories[1],")","(?=.*",venn_categories[2],")","(?=.*",venn_categories[3],")",
                                                     "(?=.*",venn_categories[4],")","(?=.*",venn_categories[5],")"),
                                    x = venn_summary$Name,perl = TRUE),][["Total"]])


  venn <- draw.quintuple.venn(
    area1 = area_n1,
    area2 = area_n2,
    area3 = area_n3,
    area4 = area_n4,
    area5 = area_n5,
    n12 = area_n12,
    n13 = area_n13,
    n14 = area_n14,
    n15 = area_n15,
    n23 = area_n23,
    n24 = area_n24,
    n25 = area_n25,
    n34 = area_n34,
    n35 = area_n35,
    n45 = area_n45,
    n123 = area_n123,
    n124 = area_n124,
    n125 = area_n125,
    n134 = area_n134,
    n135 = area_n135,
    n145 = area_n145,
    n234 = area_n234,
    n235 = area_n235,
    n245 = area_n245,
    n345 = area_n345,
    n1234 = area_n1234,
    n1235 = area_n1235,
    n1245 = area_n1245,
    n1345 = area_n1345,
    n2345 = area_n2345,
    n12345 = area_n12345,
    category = gsub(pattern = ".bed",replacement = "",x = venn_categories),
    fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
    # cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
    cat.dist = 0.25,
    cat.cex = 1.2,
    margin = 0.1,
    cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
            1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
    ind = FALSE)
  
  pdf(plot_filepath,width = 9,height = 9)
  grid.arrange(gTree(children=venn), top=paste0(SampleID," Peak Overlap")) #, bottom="subtitle")
  dev.off()
  
}
cat("\n\n\n")
