# DiffBind Volcano Plots
Stephen Kelly  
9/13/2016  



# DiffBind ChIP-Seq Differential Peaks

DiffBind is used to determine which peaks in a ChIP-Seq experiment are differential bound between sample data sets. For this report, we have subset the standard DiffBind data to plot only the single peak per gene which is closest to the gene's start site (e.g. lowest 'shortestDistance' value).

Information about the DiffBind program can be found here:

http://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf

http://bioconductor.org/packages/release/bioc/html/DiffBind.html




```r
# CODE FOR PRODUCING THE PLOTS

# function for formatting text in the report
mycat <- function(text){
    cat(gsub(pattern = "\n", replacement = "  \n", x = text))
}

# function for the new custom plots:
DiffBind_volcano_plot_top_protein_coding_promoters <- function(diff_df, signif_p_value=0.05, 
                                                               fold_change_value=1, 
                                                               plot_colors=c("grey","red","orange","blue", "green"), 
                                                               top_gene_number=10){
    # browser() # use this for debugging
    
    # this function makes the custom DiffBind volcano plots
    # only plots protein coding promoter genes from DiffBind output
    # highlights top x fold changed genes
    
    
    # ~~~~~~ PROCESS DATA ~~~~~~~~ # 
    # get the colname order for later
    original_colnames <- colnames(diff_df)
    
    suppressPackageStartupMessages(library("data.table"))
    # GET THE CLOSEST DIFFBIND PEAK PER GENE 
    # in Diffbind subset, getthe peaks closest to each gene
    # # for each gene name, find the peak with the lowest "shortestDistance" value
    # make it a data.table type
    setDT(diff_df)
    # get min value per factor, also remove duplicates
    diff_df_min <- as.data.frame(diff_df[, .SD[which.min(shortestDistance)], by=external_gene_name])
    # fix the colname order
    diff_df_min <- diff_df_min[,original_colnames]
    
    # get sample ID's from 7th, 8th colnames
    sampleID_1 <- colnames(diff_df_min)[7]
    sampleID_2 <- colnames(diff_df_min)[8]
    
    
    # subset for significant genes
    diff_df_min_signiff <- as.data.frame(subset(diff_df_min, p.value<signif_p_value & abs(Fold)>fold_change_value))
    
    # subset for protein coding genes
    diff_df_min_signiff_protein <- as.data.frame(subset(diff_df_min_signiff, gene_biotype=="protein_coding"))
    
    # subset top up/down protein coding genes
    diff_df_min_signiff_UpFold <- diff_df_min_signiff_protein[with(diff_df_min_signiff_protein, 
                                                                   order(Fold, p.value)), 
                                                              c("Fold","p.value","external_gene_name")][1:top_gene_number,]
    diff_df_min_signiff_DnFold <- diff_df_min_signiff_protein[with(diff_df_min_signiff_protein, 
                                                                   order(-Fold, p.value)), 
                                                              c("Fold","p.value","external_gene_name")][1:top_gene_number,]
    # ~~~~~~~~~~~~~~ # 
    
    
    
    
    # ~~~~~~ MAKE PLOT ~~~~~~~~ # 
    # set the multi-pane plot layout matrix
    plot_layout_matrix<-structure(c(1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 
                                    3L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 
                                    3L, 3L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 
                                    3L, 3L, 3L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 
                                    3L, 3L, 3L, 3L), 
                                  .Dim = c(11L,4L), 
                                  .Dimnames = list(NULL, c("V1", "V2", "V3", "V4")))
    
    # subset the matrix for two panels # keep this because I like this code sample
    plot_layout_matrix2 <- plot_layout_matrix[(rowSums(plot_layout_matrix < 3) > 0), , drop = FALSE]
    # x[which(x < 3, arr.ind = TRUE)[,1],]
    
    # SET PLOT LAYOUT
    layout(plot_layout_matrix2) 
    
    
    # FIRST PLOT PANEL
    # adjust first panel margins; these need to be adjusted to intended plot size
    # currently configured for 8 inch x 8 inch plota
    # plot margins # c(bottom, left, top, right) # default is c(5, 4, 4, 2) + 0.1
    par(mar=c(1,3,1,0))
    # call blank plot to fill the first panel
    plot(1,type='n',axes=FALSE,xlab="",ylab="")
    # set up the Legend in the first panel; 1st legend only title, 2nd legend only legend
    legend("center",legend = "",title=paste0("Volcano plot: ", sampleID_1, " vs. ",sampleID_2,"\n"),cex=1.3, bty='n') 
    legend("bottom",legend=c("Not Significant",
                             paste0("p < ",signif_p_value),
                             paste0("log2(Fold Change) > ",fold_change_value),
                             paste0("p < ",signif_p_value," & log2(Fold Change) > ",fold_change_value),
                             paste0("Top log2(Fold Change) genes")),
           fill=plot_colors,bty = "n",ncol=2,cex=0.9)
    
    
    
    # SECOND PLOT PANEL
    # adjust second panel margins
    par(mar=c(6,4,0,3)+ 0.1)
    
    # add base volcano plot; all points
    with(diff_df_min, plot(Fold, -log10(p.value), pch=20, xlim=c(min(Fold)-1,max(Fold)+1),col=plot_colors[1],xlab = "log2(Fold Change)"))
    
    # Add colored points for data subsets
    with(subset(diff_df_min, p.value<signif_p_value ), points(Fold, -log10(p.value), pch=20, col=plot_colors[2]))
    with(subset(diff_df_min, abs(Fold)>fold_change_value), points(Fold, -log10(p.value), pch=20, col=plot_colors[3]))
    with(subset(diff_df_min, p.value<signif_p_value & abs(Fold)>fold_change_value), points(Fold, -log10(p.value), pch=20, col=plot_colors[4]))
    
    # add points and gene labels for top genes
    suppressPackageStartupMessages(library("calibrate"))
    # up genes
    with(na.omit(diff_df_min_signiff_UpFold), points(Fold, -log10(p.value), pch=20, col="green"))
    with(na.omit(diff_df_min_signiff_UpFold), textxy(Fold, -log10(p.value), labs=external_gene_name, cex=.8))
    # down genes
    with(na.omit(diff_df_min_signiff_DnFold), points(Fold, -log10(p.value), pch=20, col="green"))
    with(na.omit(diff_df_min_signiff_DnFold), textxy(Fold, -log10(p.value), labs=external_gene_name, cex=.8))
    
    
    
    # print some gene table stats into the report; this could use optimization
    mycat('\n\n')
    mycat('\n******\n')
    # FULL DATASET STATS
    mycat('### Overall DiffBind Stats\n')
    mycat(paste0('Total number of DiffBind peaks:\n', length(as.data.frame(diff_df)[["external_gene_name"]]), '\n\n'))
    mycat(paste0('Total number of DiffBind genes:\n', length(unique(diff_df_min[["external_gene_name"]])), '\n\n'))
    
    mycat(paste0('Total number positive fold change genes:\n', 
                 length(unique(subset(diff_df_min, Fold > 0 )[["external_gene_name"]])), '\n\n'))
    
    mycat(paste0('Total number negative fold change genes:\n', 
                 length(unique(subset(diff_df_min, Fold < 0 )[["external_gene_name"]])), '\n\n'))
    
    # significant 
    mycat('\n******\n')
    mycat(paste0('Total number of p <' ,signif_p_value,' genes:\n', 
                 length(unique(subset(diff_df_min, p.value<signif_p_value )[["external_gene_name"]])), '\n\n'))
    
    mycat(paste0('Total number of p <' ,signif_p_value,' genes (pos. FC):\n', 
                 length(unique(subset(diff_df_min, p.value<signif_p_value &
                                          Fold > 0)[["external_gene_name"]])), '\n\n'))
    
    mycat(paste0('Total number of p <' ,signif_p_value,' genes (neg. FC):\n', 
                 length(unique(subset(diff_df_min, p.value<signif_p_value &
                                          Fold < 0)[["external_gene_name"]])), '\n\n'))
    
    # fold change
    mycat('\n******\n')
    mycat(paste0("Total number of log2(Fold Change) > ",fold_change_value, ' genes:\n',
                 length(unique(subset(diff_df_min, abs(Fold)>fold_change_value)[["external_gene_name"]])), '\n\n'))
    
    mycat(paste0("Total number of log2(Fold Change) > ",fold_change_value, ' genes (pos. FC):\n',
                 length(unique(subset(diff_df_min, abs(Fold)>fold_change_value & 
                                          Fold > 0)[["external_gene_name"]])), '\n\n'))
    
    mycat(paste0("Total number of log2(Fold Change) > ",fold_change_value, ' genes (neg. FC):\n',
                 length(unique(subset(diff_df_min, abs(Fold)>fold_change_value & 
                                          Fold < 0)[["external_gene_name"]])), '\n\n'))
    
    # siginificant & fold change
    mycat('\n******\n')
    mycat(paste0("Total number of p < ",signif_p_value," & log2(Fold Change) > ",fold_change_value, ' genes:\n',
                 length(unique(subset(diff_df_min, p.value<signif_p_value &
                                          abs(Fold)>fold_change_value)[["external_gene_name"]])), '\n\n'))
    
    mycat(paste0("Total number of p < ",signif_p_value," & log2(Fold Change) > ",fold_change_value, ' genes (pos. FC):\n',
                 length(unique(subset(diff_df_min, p.value<signif_p_value & 
                                          Fold > 0 &
                                          abs(Fold)>fold_change_value)[["external_gene_name"]])), '\n\n'))
    
    mycat(paste0("Total number of p < ",signif_p_value," & log2(Fold Change) > ",fold_change_value, ' genes (neg. FC):\n',
                 length(unique(subset(diff_df_min, p.value<signif_p_value & 
                                          Fold < 0 &
                                          abs(Fold)>fold_change_value)[["external_gene_name"]])), '\n\n'))
    
    
    mycat('\n\n')
    mycat('\n******\n')
    # ONLY PROTEIN CODING GENES
    mycat('### Protein Coding Gene Stats\n') # gene_biotype=="protein_coding"
    mycat(paste0('Total number of DiffBind peaks:\n', 
                 length(subset(as.data.frame(diff_df), gene_biotype=="protein_coding")[["external_gene_name"]]), '\n\n'))
    mycat(paste0('Total number of DiffBind genes:\n', 
                 length(unique(subset(diff_df_min, gene_biotype=="protein_coding")[["external_gene_name"]])), '\n\n'))
    
    mycat(paste0('Total number positive fold change genes:\n', 
                 length(unique(subset(diff_df_min, gene_biotype=="protein_coding" & 
                                          Fold > 0 )[["external_gene_name"]])), '\n\n'))
    
    mycat(paste0('Total number negative fold change genes:\n', 
                 length(unique(subset(diff_df_min, gene_biotype=="protein_coding" &
                                          Fold < 0 )[["external_gene_name"]])), '\n\n'))
    
    # significant 
    mycat('\n******\n')
    mycat(paste0('Total number of p <' ,signif_p_value,' genes:\n', 
                 length(unique(subset(diff_df_min, gene_biotype=="protein_coding" &
                                          p.value<signif_p_value )[["external_gene_name"]])), '\n\n'))
    
    mycat(paste0('Total number of p <' ,signif_p_value,' genes (pos. FC):\n', 
                 length(unique(subset(diff_df_min, gene_biotype=="protein_coding" &
                                          p.value<signif_p_value &
                                          Fold > 0)[["external_gene_name"]])), '\n\n'))
    
    mycat(paste0('Total number of p <' ,signif_p_value,' genes (neg. FC):\n', 
                 length(unique(subset(diff_df_min, gene_biotype=="protein_coding" &
                                          p.value<signif_p_value &
                                          Fold < 0)[["external_gene_name"]])), '\n\n'))
    
    # fold change
    mycat('\n******\n')
    mycat(paste0("Total number of log2(Fold Change) > ",fold_change_value, ' genes:\n',
                 length(unique(subset(diff_df_min, gene_biotype=="protein_coding" &
                                          abs(Fold)>fold_change_value)[["external_gene_name"]])), '\n\n'))
    
    mycat(paste0("Total number of log2(Fold Change) > ",fold_change_value, ' genes (pos. FC):\n',
                 length(unique(subset(diff_df_min, gene_biotype=="protein_coding" &
                                          abs(Fold)>fold_change_value & 
                                          Fold > 0)[["external_gene_name"]])), '\n\n'))
    
    mycat(paste0("Total number of log2(Fold Change) > ",fold_change_value, ' genes (neg. FC):\n',
                 length(unique(subset(diff_df_min, gene_biotype=="protein_coding" &
                                          abs(Fold)>fold_change_value & 
                                          Fold < 0)[["external_gene_name"]])), '\n\n'))
    
    # siginificant & fold change
    mycat('\n******\n')
    mycat(paste0("Total number of p < ",signif_p_value," & log2(Fold Change) > ",fold_change_value, ' genes:\n',
                 length(unique(subset(diff_df_min, gene_biotype=="protein_coding" &
                                          p.value<signif_p_value &
                                          abs(Fold)>fold_change_value)[["external_gene_name"]])), '\n\n'))
    
    mycat(paste0("Total number of p < ",signif_p_value," & log2(Fold Change) > ",fold_change_value, ' genes (pos. FC):\n',
                 length(unique(subset(diff_df_min, gene_biotype=="protein_coding" &
                                          p.value<signif_p_value & 
                                          Fold > 0 &
                                          abs(Fold)>fold_change_value)[["external_gene_name"]])), '\n\n'))
    
    mycat(paste0("Total number of p < ",signif_p_value," & log2(Fold Change) > ",fold_change_value, ' genes (neg. FC):\n',
                 length(unique(subset(diff_df_min, gene_biotype=="protein_coding" &
                                          p.value<signif_p_value & 
                                          Fold < 0 &
                                          abs(Fold)>fold_change_value)[["external_gene_name"]])), '\n\n'))
    
    
    
    mycat('\n')
    mycat('\n******\n')
    mycat('\n******\n')
    
    
}
```


```r
mycat('# DiffBind Plots and Results {.tabset}\n\n') # .tabset-fade .tabset-pills
```

# DiffBind Plots and Results {.tabset}  
  

```r
out_dir <- "/ifs/home/kellys04/projects/Bioinformatics/DiffBind_scripts_reports/DiffBind_Volcano_Plot_report/input"

mycat(paste0("Project dir:\n",out_dir,'\n'))
```

Project dir:  
/ifs/home/kellys04/projects/Bioinformatics/DiffBind_scripts_reports/DiffBind_Volcano_Plot_report/input  

```r
mycat('\n******\n')
```

  
******  

```r
diff_cl4_ctrl_file <- paste0(out_dir,"/diff_bind.Treatment4-ChIPSeq-vs-Control-ChIPSeq.p100.csv")
diff_cl4_cl5_file <- paste0(out_dir,"/diff_bind.Treatment4-ChIPSeq-vs-Treatment5-ChIPSeq.p100.csv")
diff_cl5_ctrl_file <- paste0(out_dir,"/diff_bind.Treatment5-ChIPSeq-vs-Control-ChIPSeq.p100.csv")

# pdf(file = paste0(diff_cl4_ctrl_file,"_volcano_plot.pdf"),height = 8,width = 8)
mycat(paste0("## ", basename(diff_cl4_ctrl_file), ' {.tabset}\n'))
```

## diff_bind.Treatment4-ChIPSeq-vs-Control-ChIPSeq.p100.csv {.tabset}  

```r
diff_df <- read.delim(file = diff_cl4_ctrl_file,header = TRUE,sep = ',')
DiffBind_volcano_plot_top_protein_coding_promoters(diff_df = diff_df)
```

![](DiffBind_Volcano_Plots_files/figure-html/unnamed-chunk-2-1.png)<!-- -->  
  
  
******  
### Overall DiffBind Stats  
Total number of DiffBind peaks:  
3192  
  
Total number of DiffBind genes:  
1824  
  
Total number positive fold change genes:  
1474  
  
Total number negative fold change genes:  
350  
  
  
******  
Total number of p <0.05 genes:  
1712  
  
Total number of p <0.05 genes (pos. FC):  
1377  
  
Total number of p <0.05 genes (neg. FC):  
335  
  
  
******  
Total number of log2(Fold Change) > 1 genes:  
1691  
  
Total number of log2(Fold Change) > 1 genes (pos. FC):  
1361  
  
Total number of log2(Fold Change) > 1 genes (neg. FC):  
330  
  
  
******  
Total number of p < 0.05 & log2(Fold Change) > 1 genes:  
1679  
  
Total number of p < 0.05 & log2(Fold Change) > 1 genes (pos. FC):  
1349  
  
Total number of p < 0.05 & log2(Fold Change) > 1 genes (neg. FC):  
330  
  
  
  
  
******  
### Protein Coding Gene Stats  
Total number of DiffBind peaks:  
2030  
  
Total number of DiffBind genes:  
1113  
  
Total number positive fold change genes:  
918  
  
Total number negative fold change genes:  
195  
  
  
******  
Total number of p <0.05 genes:  
1037  
  
Total number of p <0.05 genes (pos. FC):  
851  
  
Total number of p <0.05 genes (neg. FC):  
186  
  
  
******  
Total number of log2(Fold Change) > 1 genes:  
1021  
  
Total number of log2(Fold Change) > 1 genes (pos. FC):  
840  
  
Total number of log2(Fold Change) > 1 genes (neg. FC):  
181  
  
  
******  
Total number of p < 0.05 & log2(Fold Change) > 1 genes:  
1011  
  
Total number of p < 0.05 & log2(Fold Change) > 1 genes (pos. FC):  
830  
  
Total number of p < 0.05 & log2(Fold Change) > 1 genes (neg. FC):  
181  
  
  
  
******  
  
******  

```r
# dev.off()

# pdf(file = paste0(diff_cl4_cl5_file,"_volcano_plot.pdf"),height = 8,width = 8)
mycat(paste0("## ", basename(diff_cl4_cl5_file), ' {.tabset}\n'))
```

## diff_bind.Treatment4-ChIPSeq-vs-Treatment5-ChIPSeq.p100.csv {.tabset}  

```r
diff_df <- read.delim(file = diff_cl4_cl5_file,header = TRUE,sep = ',')
DiffBind_volcano_plot_top_protein_coding_promoters(diff_df = diff_df)
```

![](DiffBind_Volcano_Plots_files/figure-html/unnamed-chunk-2-2.png)<!-- -->  
  
  
******  
### Overall DiffBind Stats  
Total number of DiffBind peaks:  
279  
  
Total number of DiffBind genes:  
207  
  
Total number positive fold change genes:  
19  
  
Total number negative fold change genes:  
187  
  
  
******  
Total number of p <0.05 genes:  
37  
  
Total number of p <0.05 genes (pos. FC):  
0  
  
Total number of p <0.05 genes (neg. FC):  
37  
  
  
******  
Total number of log2(Fold Change) > 1 genes:  
0  
  
Total number of log2(Fold Change) > 1 genes (pos. FC):  
0  
  
Total number of log2(Fold Change) > 1 genes (neg. FC):  
0  
  
  
******  
Total number of p < 0.05 & log2(Fold Change) > 1 genes:  
0  
  
Total number of p < 0.05 & log2(Fold Change) > 1 genes (pos. FC):  
0  
  
Total number of p < 0.05 & log2(Fold Change) > 1 genes (neg. FC):  
0  
  
  
  
  
******  
### Protein Coding Gene Stats  
Total number of DiffBind peaks:  
176  
  
Total number of DiffBind genes:  
126  
  
Total number positive fold change genes:  
9  
  
Total number negative fold change genes:  
116  
  
  
******  
Total number of p <0.05 genes:  
22  
  
Total number of p <0.05 genes (pos. FC):  
0  
  
Total number of p <0.05 genes (neg. FC):  
22  
  
  
******  
Total number of log2(Fold Change) > 1 genes:  
0  
  
Total number of log2(Fold Change) > 1 genes (pos. FC):  
0  
  
Total number of log2(Fold Change) > 1 genes (neg. FC):  
0  
  
  
******  
Total number of p < 0.05 & log2(Fold Change) > 1 genes:  
0  
  
Total number of p < 0.05 & log2(Fold Change) > 1 genes (pos. FC):  
0  
  
Total number of p < 0.05 & log2(Fold Change) > 1 genes (neg. FC):  
0  
  
  
  
******  
  
******  

```r
# dev.off()

# pdf(file = paste0(diff_cl5_ctrl_file,"_volcano_plot.pdf"),height = 8,width = 8)
mycat(paste0("## ", basename(diff_cl5_ctrl_file), ' {.tabset}\n'))
```

## diff_bind.Treatment5-ChIPSeq-vs-Control-ChIPSeq.p100.csv {.tabset}  

```r
diff_df <- read.delim(file = diff_cl5_ctrl_file,header = TRUE,sep = ',')
DiffBind_volcano_plot_top_protein_coding_promoters(diff_df = diff_df)
```

![](DiffBind_Volcano_Plots_files/figure-html/unnamed-chunk-2-3.png)<!-- -->  
  
  
******  
### Overall DiffBind Stats  
Total number of DiffBind peaks:  
3192  
  
Total number of DiffBind genes:  
1824  
  
Total number positive fold change genes:  
1485  
  
Total number negative fold change genes:  
339  
  
  
******  
Total number of p <0.05 genes:  
1782  
  
Total number of p <0.05 genes (pos. FC):  
1452  
  
Total number of p <0.05 genes (neg. FC):  
330  
  
  
******  
Total number of log2(Fold Change) > 1 genes:  
1751  
  
Total number of log2(Fold Change) > 1 genes (pos. FC):  
1425  
  
Total number of log2(Fold Change) > 1 genes (neg. FC):  
326  
  
  
******  
Total number of p < 0.05 & log2(Fold Change) > 1 genes:  
1751  
  
Total number of p < 0.05 & log2(Fold Change) > 1 genes (pos. FC):  
1425  
  
Total number of p < 0.05 & log2(Fold Change) > 1 genes (neg. FC):  
326  
  
  
  
  
******  
### Protein Coding Gene Stats  
Total number of DiffBind peaks:  
2030  
  
Total number of DiffBind genes:  
1113  
  
Total number positive fold change genes:  
926  
  
Total number negative fold change genes:  
187  
  
  
******  
Total number of p <0.05 genes:  
1087  
  
Total number of p <0.05 genes (pos. FC):  
906  
  
Total number of p <0.05 genes (neg. FC):  
181  
  
  
******  
Total number of log2(Fold Change) > 1 genes:  
1066  
  
Total number of log2(Fold Change) > 1 genes (pos. FC):  
889  
  
Total number of log2(Fold Change) > 1 genes (neg. FC):  
177  
  
  
******  
Total number of p < 0.05 & log2(Fold Change) > 1 genes:  
1066  
  
Total number of p < 0.05 & log2(Fold Change) > 1 genes (pos. FC):  
889  
  
Total number of p < 0.05 & log2(Fold Change) > 1 genes (neg. FC):  
177  
  
  
  
******  
  
******  

```r
# dev.off()
```


## System Information


```r
system('uname -srv',intern=T)
```

```
## [1] "Linux 2.6.32-642.3.1.el6.x86_64 #1 SMP Tue Jul 12 18:30:56 UTC 2016"
```

```r
sessionInfo()
```

```
## R version 3.3.1 (2016-06-21)
## Platform: x86_64-redhat-linux-gnu (64-bit)
## Running under: CentOS release 6.8 (Final)
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] calibrate_1.7.2  MASS_7.3-45      data.table_1.9.6
## 
## loaded via a namespace (and not attached):
##  [1] magrittr_1.5    formatR_1.4     tools_3.3.1     htmltools_0.3.5
##  [5] yaml_2.1.13     Rcpp_0.12.6     stringi_1.1.1   rmarkdown_1.0  
##  [9] knitr_1.13      stringr_1.0.0   digest_0.6.10   chron_2.3-47   
## [13] evaluate_0.9
```

```r
# save.image(compress = TRUE, )
```
