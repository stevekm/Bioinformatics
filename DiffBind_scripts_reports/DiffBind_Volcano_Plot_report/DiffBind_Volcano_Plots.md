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

# Plotly Test


```r
# https://plot.ly/r/knitr/
# install.packages("plotly")
suppressPackageStartupMessages(library("plotly"))
# p <- plot_ly(economics, x = date, y = unemploy / pop)
# p
p <- plot_ly(economics, x = date, y = uempmed, name = "unemployment")
p <- p %>% add_trace(y = fitted(loess(uempmed ~ as.numeric(date))), x = date)

p2 <- plot_ly(data = iris, x = Sepal.Length, y = Petal.Length, mode = "markers", color = Species)

htmltools::tagList(list(as.widget(p), as.widget(p2)))
```

<!--html_preserve--><div id="htmlwidget-7bb77ebfde0e2e5c66a0" style="width:100%;height:400px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-7bb77ebfde0e2e5c66a0">{"x":{"data":[{"type":"scatter","inherit":false,"x":["1967-07-01","1967-08-01","1967-09-01","1967-10-01","1967-11-01","1967-12-01","1968-01-01","1968-02-01","1968-03-01","1968-04-01","1968-05-01","1968-06-01","1968-07-01","1968-08-01","1968-09-01","1968-10-01","1968-11-01","1968-12-01","1969-01-01","1969-02-01","1969-03-01","1969-04-01","1969-05-01","1969-06-01","1969-07-01","1969-08-01","1969-09-01","1969-10-01","1969-11-01","1969-12-01","1970-01-01","1970-02-01","1970-03-01","1970-04-01","1970-05-01","1970-06-01","1970-07-01","1970-08-01","1970-09-01","1970-10-01","1970-11-01","1970-12-01","1971-01-01","1971-02-01","1971-03-01","1971-04-01","1971-05-01","1971-06-01","1971-07-01","1971-08-01","1971-09-01","1971-10-01","1971-11-01","1971-12-01","1972-01-01","1972-02-01","1972-03-01","1972-04-01","1972-05-01","1972-06-01","1972-07-01","1972-08-01","1972-09-01","1972-10-01","1972-11-01","1972-12-01","1973-01-01","1973-02-01","1973-03-01","1973-04-01","1973-05-01","1973-06-01","1973-07-01","1973-08-01","1973-09-01","1973-10-01","1973-11-01","1973-12-01","1974-01-01","1974-02-01","1974-03-01","1974-04-01","1974-05-01","1974-06-01","1974-07-01","1974-08-01","1974-09-01","1974-10-01","1974-11-01","1974-12-01","1975-01-01","1975-02-01","1975-03-01","1975-04-01","1975-05-01","1975-06-01","1975-07-01","1975-08-01","1975-09-01","1975-10-01","1975-11-01","1975-12-01","1976-01-01","1976-02-01","1976-03-01","1976-04-01","1976-05-01","1976-06-01","1976-07-01","1976-08-01","1976-09-01","1976-10-01","1976-11-01","1976-12-01","1977-01-01","1977-02-01","1977-03-01","1977-04-01","1977-05-01","1977-06-01","1977-07-01","1977-08-01","1977-09-01","1977-10-01","1977-11-01","1977-12-01","1978-01-01","1978-02-01","1978-03-01","1978-04-01","1978-05-01","1978-06-01","1978-07-01","1978-08-01","1978-09-01","1978-10-01","1978-11-01","1978-12-01","1979-01-01","1979-02-01","1979-03-01","1979-04-01","1979-05-01","1979-06-01","1979-07-01","1979-08-01","1979-09-01","1979-10-01","1979-11-01","1979-12-01","1980-01-01","1980-02-01","1980-03-01","1980-04-01","1980-05-01","1980-06-01","1980-07-01","1980-08-01","1980-09-01","1980-10-01","1980-11-01","1980-12-01","1981-01-01","1981-02-01","1981-03-01","1981-04-01","1981-05-01","1981-06-01","1981-07-01","1981-08-01","1981-09-01","1981-10-01","1981-11-01","1981-12-01","1982-01-01","1982-02-01","1982-03-01","1982-04-01","1982-05-01","1982-06-01","1982-07-01","1982-08-01","1982-09-01","1982-10-01","1982-11-01","1982-12-01","1983-01-01","1983-02-01","1983-03-01","1983-04-01","1983-05-01","1983-06-01","1983-07-01","1983-08-01","1983-09-01","1983-10-01","1983-11-01","1983-12-01","1984-01-01","1984-02-01","1984-03-01","1984-04-01","1984-05-01","1984-06-01","1984-07-01","1984-08-01","1984-09-01","1984-10-01","1984-11-01","1984-12-01","1985-01-01","1985-02-01","1985-03-01","1985-04-01","1985-05-01","1985-06-01","1985-07-01","1985-08-01","1985-09-01","1985-10-01","1985-11-01","1985-12-01","1986-01-01","1986-02-01","1986-03-01","1986-04-01","1986-05-01","1986-06-01","1986-07-01","1986-08-01","1986-09-01","1986-10-01","1986-11-01","1986-12-01","1987-01-01","1987-02-01","1987-03-01","1987-04-01","1987-05-01","1987-06-01","1987-07-01","1987-08-01","1987-09-01","1987-10-01","1987-11-01","1987-12-01","1988-01-01","1988-02-01","1988-03-01","1988-04-01","1988-05-01","1988-06-01","1988-07-01","1988-08-01","1988-09-01","1988-10-01","1988-11-01","1988-12-01","1989-01-01","1989-02-01","1989-03-01","1989-04-01","1989-05-01","1989-06-01","1989-07-01","1989-08-01","1989-09-01","1989-10-01","1989-11-01","1989-12-01","1990-01-01","1990-02-01","1990-03-01","1990-04-01","1990-05-01","1990-06-01","1990-07-01","1990-08-01","1990-09-01","1990-10-01","1990-11-01","1990-12-01","1991-01-01","1991-02-01","1991-03-01","1991-04-01","1991-05-01","1991-06-01","1991-07-01","1991-08-01","1991-09-01","1991-10-01","1991-11-01","1991-12-01","1992-01-01","1992-02-01","1992-03-01","1992-04-01","1992-05-01","1992-06-01","1992-07-01","1992-08-01","1992-09-01","1992-10-01","1992-11-01","1992-12-01","1993-01-01","1993-02-01","1993-03-01","1993-04-01","1993-05-01","1993-06-01","1993-07-01","1993-08-01","1993-09-01","1993-10-01","1993-11-01","1993-12-01","1994-01-01","1994-02-01","1994-03-01","1994-04-01","1994-05-01","1994-06-01","1994-07-01","1994-08-01","1994-09-01","1994-10-01","1994-11-01","1994-12-01","1995-01-01","1995-02-01","1995-03-01","1995-04-01","1995-05-01","1995-06-01","1995-07-01","1995-08-01","1995-09-01","1995-10-01","1995-11-01","1995-12-01","1996-01-01","1996-02-01","1996-03-01","1996-04-01","1996-05-01","1996-06-01","1996-07-01","1996-08-01","1996-09-01","1996-10-01","1996-11-01","1996-12-01","1997-01-01","1997-02-01","1997-03-01","1997-04-01","1997-05-01","1997-06-01","1997-07-01","1997-08-01","1997-09-01","1997-10-01","1997-11-01","1997-12-01","1998-01-01","1998-02-01","1998-03-01","1998-04-01","1998-05-01","1998-06-01","1998-07-01","1998-08-01","1998-09-01","1998-10-01","1998-11-01","1998-12-01","1999-01-01","1999-02-01","1999-03-01","1999-04-01","1999-05-01","1999-06-01","1999-07-01","1999-08-01","1999-09-01","1999-10-01","1999-11-01","1999-12-01","2000-01-01","2000-02-01","2000-03-01","2000-04-01","2000-05-01","2000-06-01","2000-07-01","2000-08-01","2000-09-01","2000-10-01","2000-11-01","2000-12-01","2001-01-01","2001-02-01","2001-03-01","2001-04-01","2001-05-01","2001-06-01","2001-07-01","2001-08-01","2001-09-01","2001-10-01","2001-11-01","2001-12-01","2002-01-01","2002-02-01","2002-03-01","2002-04-01","2002-05-01","2002-06-01","2002-07-01","2002-08-01","2002-09-01","2002-10-01","2002-11-01","2002-12-01","2003-01-01","2003-02-01","2003-03-01","2003-04-01","2003-05-01","2003-06-01","2003-07-01","2003-08-01","2003-09-01","2003-10-01","2003-11-01","2003-12-01","2004-01-01","2004-02-01","2004-03-01","2004-04-01","2004-05-01","2004-06-01","2004-07-01","2004-08-01","2004-09-01","2004-10-01","2004-11-01","2004-12-01","2005-01-01","2005-02-01","2005-03-01","2005-04-01","2005-05-01","2005-06-01","2005-07-01","2005-08-01","2005-09-01","2005-10-01","2005-11-01","2005-12-01","2006-01-01","2006-02-01","2006-03-01","2006-04-01","2006-05-01","2006-06-01","2006-07-01","2006-08-01","2006-09-01","2006-10-01","2006-11-01","2006-12-01","2007-01-01","2007-02-01","2007-03-01","2007-04-01","2007-05-01","2007-06-01","2007-07-01","2007-08-01","2007-09-01","2007-10-01","2007-11-01","2007-12-01","2008-01-01","2008-02-01","2008-03-01","2008-04-01","2008-05-01","2008-06-01","2008-07-01","2008-08-01","2008-09-01","2008-10-01","2008-11-01","2008-12-01","2009-01-01","2009-02-01","2009-03-01","2009-04-01","2009-05-01","2009-06-01","2009-07-01","2009-08-01","2009-09-01","2009-10-01","2009-11-01","2009-12-01","2010-01-01","2010-02-01","2010-03-01","2010-04-01","2010-05-01","2010-06-01","2010-07-01","2010-08-01","2010-09-01","2010-10-01","2010-11-01","2010-12-01","2011-01-01","2011-02-01","2011-03-01","2011-04-01","2011-05-01","2011-06-01","2011-07-01","2011-08-01","2011-09-01","2011-10-01","2011-11-01","2011-12-01","2012-01-01","2012-02-01","2012-03-01","2012-04-01","2012-05-01","2012-06-01","2012-07-01","2012-08-01","2012-09-01","2012-10-01","2012-11-01","2012-12-01","2013-01-01","2013-02-01","2013-03-01","2013-04-01","2013-05-01","2013-06-01","2013-07-01","2013-08-01","2013-09-01","2013-10-01","2013-11-01","2013-12-01","2014-01-01","2014-02-01","2014-03-01","2014-04-01","2014-05-01","2014-06-01","2014-07-01","2014-08-01","2014-09-01","2014-10-01","2014-11-01","2014-12-01","2015-01-01","2015-02-01","2015-03-01","2015-04-01"],"y":[4.5,4.7,4.6,4.9,4.7,4.8,5.1,4.5,4.1,4.6,4.4,4.4,4.5,4.2,4.6,4.8,4.4,4.4,4.4,4.9,4,4,4.2,4.4,4.4,4.4,4.7,4.5,4.8,4.6,4.6,4.5,4.6,4.1,4.7,4.9,5.1,5.4,5.2,5.2,5.6,5.9,6.2,6.3,6.4,6.5,6.7,5.7,6.2,6.4,5.8,6.5,6.4,6.2,6.2,6.6,6.6,6.7,6.6,5.4,6.1,6,5.6,5.7,5.7,6.1,5.7,5.2,5.5,5,4.9,5,5.2,4.9,5.4,5.5,5.1,4.7,5,5.1,4.8,5,4.6,5.3,5.7,5,5.3,5.5,5.2,5.7,6.3,7.1,7.2,8.7,9.4,8.8,8.6,9.2,9.2,8.6,9.5,9,9,8.2,8.7,8.2,8.3,7.8,7.7,7.9,7.8,7.7,8.4,8,7.5,7.2,7.2,7.3,7.9,6.2,7.1,7,6.7,6.9,7,6.8,6.5,6.7,6.2,6.1,5.7,6,5.8,5.8,5.6,5.9,5.5,5.6,5.9,5.9,5.9,5.4,5.6,5.6,5.9,4.8,5.5,5.5,5.3,5.7,5.3,5.8,6,5.8,5.7,6.4,7,7.5,7.7,7.5,7.7,7.5,7.4,7.1,7.1,7.4,6.9,6.6,7.1,7.2,6.8,6.8,6.9,6.9,7.1,7.5,7.7,8.1,8.5,9.5,8.5,8.7,9.5,9.7,10,10.2,11.1,9.8,10.4,10.9,12.3,11.3,10.1,9.3,9.3,9.4,9.3,8.7,9.1,8.3,8.3,8.2,9.1,7.5,7.5,7.3,7.6,7.2,7.2,7.3,6.8,7.1,7.1,6.9,6.9,6.6,6.9,7.1,6.9,7.1,7,6.8,6.7,6.9,6.8,6.7,6.8,7,6.9,7.1,7.4,7,7.1,7.1,6.9,6.6,6.6,7.1,6.6,6.5,6.5,6.4,6,6.3,6.2,6,6.2,6.3,6.4,5.9,5.9,5.8,6.1,5.9,5.7,5.6,5.7,5.9,5.6,5.4,5.4,5.4,5.3,5.4,5.6,5,4.9,4.9,4.8,4.9,5.1,5.3,5.1,4.8,5.2,5.2,5.4,5.4,5.6,5.8,5.7,5.9,6,6.2,6.7,6.6,6.4,6.9,7,7.3,6.8,7.2,7.5,7.8,8.1,8.2,8.3,8.5,8.8,8.7,8.6,8.8,8.6,9,9,9.3,8.6,8.5,8.5,8.4,8.1,8.3,8.2,8.2,8.3,8,8.3,8.3,8.6,9.2,9.3,9.1,9.2,9.3,9,8.9,9.2,10,9,8.7,8,8.1,8.3,8.3,9.1,7.9,8.5,8.3,7.9,8.2,8,8.3,8.3,7.8,8.3,8.6,8.6,8.3,8.3,8.4,8.5,8.3,7.7,7.8,7.8,8.1,7.9,8.3,8,8,8.3,7.8,8.2,7.7,7.6,7.5,7.4,7,6.8,6.7,6,6.9,6.7,6.8,6.7,5.8,6.6,6.8,6.9,6.8,6.8,6.2,6.5,6.3,5.8,6.5,6,6.1,6.2,5.8,5.8,6.1,6,6.1,5.8,5.7,6,6.3,5.2,6.1,6.1,6,5.8,6.1,6.6,5.9,6.3,6,6.8,6.9,7.2,7.3,7.7,8.2,8.4,8.3,8.4,8.9,9.5,11,8.9,9,9.5,9.6,9.3,9.6,9.6,9.5,9.7,10.2,9.9,11.5,10.3,10.1,10.2,10.4,10.3,10.4,10.6,10.2,10.2,9.5,9.9,11,8.9,9.2,9.6,9.5,9.7,9.5,9.4,9.2,9.3,9,9.1,9,8.8,9.2,8.4,8.6,8.5,8.7,8.6,9.1,8.7,8.4,8.5,7.3,8,8.4,8,7.9,8.3,7.5,8.3,8.5,9.1,8.6,8.2,7.7,8.7,8.8,8.7,8.4,8.6,8.4,9,8.7,8.7,9.4,7.9,9,9.7,9.7,10.2,10.4,9.8,10.5,10.7,11.7,12.3,13.1,14.2,17.2,16,16.3,17.8,18.9,19.8,20.1,20,19.9,20.4,22.1,22.3,25.2,22.3,21,20.3,21.2,21,21.9,21.6,21.1,21.5,20.9,21.6,22.3,22,22.4,22,20.5,20.9,20.5,21,19.8,19.2,19.1,19.9,20.1,17.5,18.5,18.8,19.7,18.5,17.6,16.2,17.5,17.7,17.1,17,16.6,16.3,16.8,16.5,16.1,17,17,15.9,16.2,15.9,15.6,14.5,13.2,13.5,13.3,13.3,13.5,12.8,12.6,13.4,13.1,12.2,11.7],"name":"unemployment"},{"y":[4.43762517079951,4.4645618436822,4.49138942311142,4.51724708919133,4.54385802906226,4.56950477874467,4.59589646428572,4.62217572343544,4.64665716790769,4.67271690920106,4.6978270882456,4.72366105971261,4.74855148661636,4.77415707348543,4.79964559399087,4.82419980532375,4.84945599764127,4.87378411003184,4.89880535938984,4.92370620944637,4.94609325234004,4.97076308965266,4.99452062979751,5.01894905151925,5.04247169846085,5.06665608982011,5.0907154636895,5.1138790593595,5.13769045927644,5.16061281494168,5.18417362613457,5.20760608689975,5.22866000852895,5.25184632376158,5.27416062684618,5.2970898807234,5.31915414954231,5.34182372729251,5.36436033642662,5.38604271827318,5.4083157076903,5.42974170846965,5.45174846339819,5.47361891677278,5.49325515885432,5.51486433390772,5.53564480177152,5.55698126970497,5.57749656224063,5.59855770828254,5.61947793458211,5.63958850444468,5.66022946526274,5.68006851299565,5.70042759356059,5.72064242144536,5.73942187572145,5.75935558648111,5.77850704601405,5.79815236102412,5.81702346268411,5.836377767173,5.85558317900941,5.87402665973199,5.89293711745726,5.91109389350537,5.92970678234788,5.94816744560006,5.96471019019949,5.98287926809349,6.00031544787037,6.01818048425026,6.03531294181054,6.05284707752503,6.07021071128564,6.08685348581497,6.10388657125709,6.12021264005503,6.1369217486374,6.15346873162407,6.16827645496633,6.18451956518863,6.20008901299332,6.21602441340947,6.2312992572689,6.24693351933765,6.2624172622678,6.27725981390812,6.2924523295434,6.30701622805582,6.32192408785497,6.3366898048737,6.34990581848157,6.36440569454565,6.37830716785381,6.39253865513537,6.40618357710178,6.42015324705955,6.433992380694,6.44726240125313,6.46084962861632,6.47387904865105,6.48722094120173,6.50044067378727,6.51269843154964,6.52568654679971,6.53814410420033,6.55090325942032,6.56314242233639,6.57567918898504,6.58810545687309,6.60002726642883,6.61224105679806,6.62396042268661,6.63596830699578,6.64787406890254,6.65854129561659,6.67025738438488,6.68150297471781,6.69302942447965,6.70409467293152,6.71543805514564,6.72669092141444,6.73749620731782,6.74857591709142,6.7592168119311,6.77012993666826,6.78096092181827,6.79067455808547,6.80135390190676,6.81161521697981,6.82214424281826,6.83226326861376,6.84264854792811,6.8529632941124,6.86287974817124,6.87306065888397,6.88285077448261,6.8929044211825,6.90289591111054,6.91186846928485,6.92174634969388,6.93125108131489,6.9410179647647,6.95028774025756,6.959543790221,6.96847860738075,6.97682579958723,6.98514851273727,6.9929161788425,7.00065330872633,7.00810301066709,7.01481794357069,7.02173079093673,7.02816604715853,7.03455939573394,7.04050474380844,7.04640511433673,7.05206511859065,7.05732013272422,7.0625271574294,7.06735663687555,7.07213720277548,7.07671120726167,7.08067071946563,7.08487051396303,7.0887570204868,7.09259600588825,7.09614617493698,7.09965087198597,7.10299584834314,7.10608729823534,7.10913798473843,7.11195747061547,7.11474038800798,7.11739738956929,7.11969473185898,7.12213101018536,7.12438767349266,7.12662183241776,7.12869572892229,7.13075428918961,7.13273377434763,7.13458056648001,7.13642445147864,7.13815285026813,7.13988765585084,7.14157719118489,7.14307000960976,7.14469230846274,7.14623803503514,7.14781690418151,7.1493334346234,7.15089539480665,7.15245892546315,7.15397996631729,7.15556658650906,7.15712280469256,7.15875903516309,7.16043064096749,7.16203249418781,7.16379202595352,7.16554755427741,7.16742278079747,7.16930310440468,7.17132054642258,7.17342042544306,7.17553750880797,7.17781960182468,7.18012585393278,7.1826166816892,7.18522375130884,7.18768413348429,7.19053136064224,7.19341628274636,7.19653789435332,7.19949699120864,7.20228034558004,7.20479292235707,7.20697434370722,7.20897800655242,7.21068227027123,7.21220872127684,7.21350480399299,7.21448426321592,7.21536452306796,7.21601959469871,7.21650111838215,7.21678565602341,7.21690014563051,7.21684025892618,7.21662398634941,7.21624487856033,7.21573502207521,7.21506839529364,7.21426780150555,7.21343634736119,7.21240340411847,7.21129877574482,7.21005678918058,7.20876520054043,7.20734587274695,7.20584856992498,7.20433303635392,7.2027062052943,7.20108066153788,7.19935600433285,7.19759378140424,7.19591849457021,7.19410682278223,7.19234062447132,7.19051012191535,7.18874104000264,7.18692340877094,7.18512446790681,7.18340906488535,7.18167080833708,7.18002947179496,7.1783836016548,7.17679683118711,7.17542122690639,7.17396955646991,7.17264353342516,7.17136272506777,7.17021739214421,7.16913914795817,7.16817599542265,7.16736120667748,7.16664843130491,7.16609128276206,7.16566058664469,7.16538539148273,7.16527729482438,7.16532067310311,7.16553325825716,7.16593719146208,7.1665140404603,7.16731023068365,7.16831791384047,7.16950217230421,7.17094992547106,7.17257539786012,7.17449492312926,7.17666635063677,7.17885096567124,7.1815244400289,7.18437406631439,7.1870902315122,7.1890414085487,7.19038344388709,7.19106643656663,7.1911248376006,7.19058830838535,7.18951612942564,7.18786259011827,7.18568067557475,7.18318464057502,7.18005527450215,7.17659749821857,7.17260656129948,7.16836460926277,7.16361461330011,7.15851815165058,7.15328114727699,7.14758023454853,7.14180788000305,7.1356050279983,7.12918637772939,7.12322542580293,7.11646920671193,7.10979837442428,7.10279397478156,7.09593257025158,7.08878250185979,7.08159769071767,7.07463650471955,7.06746063279298,7.06055770013378,7.05349327922541,7.04652478298945,7.04033512850958,7.03362234258881,7.02729034415678,7.02094276797879,7.01501380081524,7.00913394618265,7.0035310717363,6.99839759397334,6.99341704903703,6.98893542812379,6.98467972450017,6.9808316684855,6.9777284086546,6.97472934209242,6.97228806737571,6.9702676008508,6.96882296091336,6.96788360622835,6.96753295466611,6.96777907499802,6.96866414324031,6.97015572393273,6.97237902378221,6.97532169417718,6.97874921868808,6.98315942938909,6.98818603790249,6.9941899194566,7.00080839257106,7.00850845678093,7.01710980108555,7.026316758135,7.03676889306001,7.04781654029176,7.0602219776398,7.07365936250531,7.08670535804221,7.10217991605711,7.11821075921032,7.13549733603649,7.15256361125926,7.17054899321049,7.18889203134793,7.20698549811489,7.22603728127949,7.24481981932461,7.26458693810362,7.28472011380933,7.30322096721675,7.32405553441379,7.34457124139572,7.36613745861168,7.38736440299286,7.40966886081501,7.43235101580558,7.45466270959868,7.47809304509753,7.50113197353764,7.52531708013201,7.54988828463531,7.57241490583421,7.59772564434336,7.62259221169884,7.64867403570804,7.6742900038316,7.70114950392811,7.72840674217525,7.75516501757848,7.78320987181502,7.81073354523407,7.83957260604726,7.86881780575158,7.89654540654235,7.92658151686064,7.95604015863958,7.98688689798146,8.01713320845448,8.0487971674072,8.08087896039957,8.11232529371279,8.14523424621809,8.17748424547314,8.21122694807863,8.24539588546426,8.27662592518292,8.3116107861035,8.34587762305295,8.38171254847824,8.41680521749059,8.45349679800785,8.49062625354697,8.52697650326726,8.56497255380043,8.60216463283139,8.64103386894587,8.6803493808227,8.71624516508837,8.75641474301468,8.79571812971861,8.83677720763076,8.87694458976642,8.91889975825165,8.96131284274096,9.0027953634523,9.04611447841678,9.08847699180762,9.13270872783454,9.17740678060602,9.21818144036373,9.26377170169922,9.3083399927416,9.35478414897802,9.40003880062014,9.4471221430545,9.49453217570042,9.54072508434968,9.58878161485598,9.63560251237825,9.68431038869504,9.73335113469861,9.77953041996728,9.82921898669442,9.87762511933585,9.92797711201142,9.97702761554656,10.0280478821207,10.0794096211804,10.1294408365368,10.1814781611575,10.2321655147746,10.2848832729064,10.337948682999,10.3861785865502,10.4399091159253,10.4922408615613,10.5466638251109,10.5996680141173,10.6547882597914,10.7102647198365,10.7642919594817,10.8204730133137,10.8751844636963,10.9320749592654,10.9893278486806,11.0413526625849,11.0992987284194,11.1557240617059,11.2143915699408,11.2715174190746,11.33091121766,11.3906759725016,11.4488672110916,11.5093655679462,11.5682690897503,11.629505896568,11.6911198391172,11.7470963397203,11.8094315158253,11.8701184114185,11.9332040381499,11.994619522067,12.0584604473749,12.1226870708245,12.1852102830152,12.2501995167037,12.3134630845854,12.3792197764629,12.4453683459573,12.5076051444595,12.5745157435744,12.6396445425807,12.7073346800225,12.7732202167405,12.8416947404595,12.9105697445941,12.9776054643998,13.047272166474,13.1150763913583,13.1855396393221,13.2564095471769,13.3207716511524,13.392418983705,13.4621457545346,13.5345993166222,13.6051031633433,13.6783570869171,13.7520170438154,13.8236870216246,13.8981445546366,13.9705854866736,14.0458387489761,14.1214956905365,14.1901775667852,14.2666014833855,14.3409431551964,14.4181582294432,14.4932646440621,14.5712690291324,14.6496738316166,14.7259302483818,14.805121576317,14.8821382531203,14.9621142596833,15.0424883295936,15.1154260109054,15.1965563675097,15.2754473115845,15.3573577751727,15.4370027678658,15.5196914916147,15.602775016867,15.6835538827646,15.7674128828051,15.8489413145511,15.9335739425567,16.0185990179992,16.0984934675879,16.1842767637489,16.2676651908116,16.3542175130179,16.43834926521,16.5256687666385,16.6133774382739,16.6986268905397,16.7870999560248,16.8730882506501,16.9623238631618,17.0519462918138,17.1332276621938,17.2235849715648,17.3113945159294,17.4025098006831,17.4910519759874,17.5829233893006,17.6751783569101,17.7648221033462,17.8578304141427,17.9482023087486,18.0419621159091,18.1361031232991,18.221460963609,18.3163261633723,18.4084928352151,18.5041049596983,18.5969935682912,18.6933507706712,18.7900859114369,18.8840599622198,18.9815373955098,19.0762289002728,19.1744467792641,19.2730402425744,19.3624148289787,19.4617217963164],"x":["1967-07-01","1967-08-01","1967-09-01","1967-10-01","1967-11-01","1967-12-01","1968-01-01","1968-02-01","1968-03-01","1968-04-01","1968-05-01","1968-06-01","1968-07-01","1968-08-01","1968-09-01","1968-10-01","1968-11-01","1968-12-01","1969-01-01","1969-02-01","1969-03-01","1969-04-01","1969-05-01","1969-06-01","1969-07-01","1969-08-01","1969-09-01","1969-10-01","1969-11-01","1969-12-01","1970-01-01","1970-02-01","1970-03-01","1970-04-01","1970-05-01","1970-06-01","1970-07-01","1970-08-01","1970-09-01","1970-10-01","1970-11-01","1970-12-01","1971-01-01","1971-02-01","1971-03-01","1971-04-01","1971-05-01","1971-06-01","1971-07-01","1971-08-01","1971-09-01","1971-10-01","1971-11-01","1971-12-01","1972-01-01","1972-02-01","1972-03-01","1972-04-01","1972-05-01","1972-06-01","1972-07-01","1972-08-01","1972-09-01","1972-10-01","1972-11-01","1972-12-01","1973-01-01","1973-02-01","1973-03-01","1973-04-01","1973-05-01","1973-06-01","1973-07-01","1973-08-01","1973-09-01","1973-10-01","1973-11-01","1973-12-01","1974-01-01","1974-02-01","1974-03-01","1974-04-01","1974-05-01","1974-06-01","1974-07-01","1974-08-01","1974-09-01","1974-10-01","1974-11-01","1974-12-01","1975-01-01","1975-02-01","1975-03-01","1975-04-01","1975-05-01","1975-06-01","1975-07-01","1975-08-01","1975-09-01","1975-10-01","1975-11-01","1975-12-01","1976-01-01","1976-02-01","1976-03-01","1976-04-01","1976-05-01","1976-06-01","1976-07-01","1976-08-01","1976-09-01","1976-10-01","1976-11-01","1976-12-01","1977-01-01","1977-02-01","1977-03-01","1977-04-01","1977-05-01","1977-06-01","1977-07-01","1977-08-01","1977-09-01","1977-10-01","1977-11-01","1977-12-01","1978-01-01","1978-02-01","1978-03-01","1978-04-01","1978-05-01","1978-06-01","1978-07-01","1978-08-01","1978-09-01","1978-10-01","1978-11-01","1978-12-01","1979-01-01","1979-02-01","1979-03-01","1979-04-01","1979-05-01","1979-06-01","1979-07-01","1979-08-01","1979-09-01","1979-10-01","1979-11-01","1979-12-01","1980-01-01","1980-02-01","1980-03-01","1980-04-01","1980-05-01","1980-06-01","1980-07-01","1980-08-01","1980-09-01","1980-10-01","1980-11-01","1980-12-01","1981-01-01","1981-02-01","1981-03-01","1981-04-01","1981-05-01","1981-06-01","1981-07-01","1981-08-01","1981-09-01","1981-10-01","1981-11-01","1981-12-01","1982-01-01","1982-02-01","1982-03-01","1982-04-01","1982-05-01","1982-06-01","1982-07-01","1982-08-01","1982-09-01","1982-10-01","1982-11-01","1982-12-01","1983-01-01","1983-02-01","1983-03-01","1983-04-01","1983-05-01","1983-06-01","1983-07-01","1983-08-01","1983-09-01","1983-10-01","1983-11-01","1983-12-01","1984-01-01","1984-02-01","1984-03-01","1984-04-01","1984-05-01","1984-06-01","1984-07-01","1984-08-01","1984-09-01","1984-10-01","1984-11-01","1984-12-01","1985-01-01","1985-02-01","1985-03-01","1985-04-01","1985-05-01","1985-06-01","1985-07-01","1985-08-01","1985-09-01","1985-10-01","1985-11-01","1985-12-01","1986-01-01","1986-02-01","1986-03-01","1986-04-01","1986-05-01","1986-06-01","1986-07-01","1986-08-01","1986-09-01","1986-10-01","1986-11-01","1986-12-01","1987-01-01","1987-02-01","1987-03-01","1987-04-01","1987-05-01","1987-06-01","1987-07-01","1987-08-01","1987-09-01","1987-10-01","1987-11-01","1987-12-01","1988-01-01","1988-02-01","1988-03-01","1988-04-01","1988-05-01","1988-06-01","1988-07-01","1988-08-01","1988-09-01","1988-10-01","1988-11-01","1988-12-01","1989-01-01","1989-02-01","1989-03-01","1989-04-01","1989-05-01","1989-06-01","1989-07-01","1989-08-01","1989-09-01","1989-10-01","1989-11-01","1989-12-01","1990-01-01","1990-02-01","1990-03-01","1990-04-01","1990-05-01","1990-06-01","1990-07-01","1990-08-01","1990-09-01","1990-10-01","1990-11-01","1990-12-01","1991-01-01","1991-02-01","1991-03-01","1991-04-01","1991-05-01","1991-06-01","1991-07-01","1991-08-01","1991-09-01","1991-10-01","1991-11-01","1991-12-01","1992-01-01","1992-02-01","1992-03-01","1992-04-01","1992-05-01","1992-06-01","1992-07-01","1992-08-01","1992-09-01","1992-10-01","1992-11-01","1992-12-01","1993-01-01","1993-02-01","1993-03-01","1993-04-01","1993-05-01","1993-06-01","1993-07-01","1993-08-01","1993-09-01","1993-10-01","1993-11-01","1993-12-01","1994-01-01","1994-02-01","1994-03-01","1994-04-01","1994-05-01","1994-06-01","1994-07-01","1994-08-01","1994-09-01","1994-10-01","1994-11-01","1994-12-01","1995-01-01","1995-02-01","1995-03-01","1995-04-01","1995-05-01","1995-06-01","1995-07-01","1995-08-01","1995-09-01","1995-10-01","1995-11-01","1995-12-01","1996-01-01","1996-02-01","1996-03-01","1996-04-01","1996-05-01","1996-06-01","1996-07-01","1996-08-01","1996-09-01","1996-10-01","1996-11-01","1996-12-01","1997-01-01","1997-02-01","1997-03-01","1997-04-01","1997-05-01","1997-06-01","1997-07-01","1997-08-01","1997-09-01","1997-10-01","1997-11-01","1997-12-01","1998-01-01","1998-02-01","1998-03-01","1998-04-01","1998-05-01","1998-06-01","1998-07-01","1998-08-01","1998-09-01","1998-10-01","1998-11-01","1998-12-01","1999-01-01","1999-02-01","1999-03-01","1999-04-01","1999-05-01","1999-06-01","1999-07-01","1999-08-01","1999-09-01","1999-10-01","1999-11-01","1999-12-01","2000-01-01","2000-02-01","2000-03-01","2000-04-01","2000-05-01","2000-06-01","2000-07-01","2000-08-01","2000-09-01","2000-10-01","2000-11-01","2000-12-01","2001-01-01","2001-02-01","2001-03-01","2001-04-01","2001-05-01","2001-06-01","2001-07-01","2001-08-01","2001-09-01","2001-10-01","2001-11-01","2001-12-01","2002-01-01","2002-02-01","2002-03-01","2002-04-01","2002-05-01","2002-06-01","2002-07-01","2002-08-01","2002-09-01","2002-10-01","2002-11-01","2002-12-01","2003-01-01","2003-02-01","2003-03-01","2003-04-01","2003-05-01","2003-06-01","2003-07-01","2003-08-01","2003-09-01","2003-10-01","2003-11-01","2003-12-01","2004-01-01","2004-02-01","2004-03-01","2004-04-01","2004-05-01","2004-06-01","2004-07-01","2004-08-01","2004-09-01","2004-10-01","2004-11-01","2004-12-01","2005-01-01","2005-02-01","2005-03-01","2005-04-01","2005-05-01","2005-06-01","2005-07-01","2005-08-01","2005-09-01","2005-10-01","2005-11-01","2005-12-01","2006-01-01","2006-02-01","2006-03-01","2006-04-01","2006-05-01","2006-06-01","2006-07-01","2006-08-01","2006-09-01","2006-10-01","2006-11-01","2006-12-01","2007-01-01","2007-02-01","2007-03-01","2007-04-01","2007-05-01","2007-06-01","2007-07-01","2007-08-01","2007-09-01","2007-10-01","2007-11-01","2007-12-01","2008-01-01","2008-02-01","2008-03-01","2008-04-01","2008-05-01","2008-06-01","2008-07-01","2008-08-01","2008-09-01","2008-10-01","2008-11-01","2008-12-01","2009-01-01","2009-02-01","2009-03-01","2009-04-01","2009-05-01","2009-06-01","2009-07-01","2009-08-01","2009-09-01","2009-10-01","2009-11-01","2009-12-01","2010-01-01","2010-02-01","2010-03-01","2010-04-01","2010-05-01","2010-06-01","2010-07-01","2010-08-01","2010-09-01","2010-10-01","2010-11-01","2010-12-01","2011-01-01","2011-02-01","2011-03-01","2011-04-01","2011-05-01","2011-06-01","2011-07-01","2011-08-01","2011-09-01","2011-10-01","2011-11-01","2011-12-01","2012-01-01","2012-02-01","2012-03-01","2012-04-01","2012-05-01","2012-06-01","2012-07-01","2012-08-01","2012-09-01","2012-10-01","2012-11-01","2012-12-01","2013-01-01","2013-02-01","2013-03-01","2013-04-01","2013-05-01","2013-06-01","2013-07-01","2013-08-01","2013-09-01","2013-10-01","2013-11-01","2013-12-01","2014-01-01","2014-02-01","2014-03-01","2014-04-01","2014-05-01","2014-06-01","2014-07-01","2014-08-01","2014-09-01","2014-10-01","2014-11-01","2014-12-01","2015-01-01","2015-02-01","2015-03-01","2015-04-01"]}],"layout":{"xaxis":{"title":"date"},"yaxis":{"title":"uempmed"},"margin":{"b":40,"l":60,"t":25,"r":10}},"url":null,"width":null,"height":null,"source":"A","config":{"modeBarButtonsToRemove":["sendDataToCloud"]},"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
<div id="htmlwidget-cb6aa5fb39eebb3b6c72" style="width:100%;height:400px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-cb6aa5fb39eebb3b6c72">{"x":{"data":[{"type":"scatter","inherit":false,"x":[6.3,5.8,7.1,6.3,6.5,7.6,4.9,7.3,6.7,7.2,6.5,6.4,6.8,5.7,5.8,6.4,6.5,7.7,7.7,6,6.9,5.6,7.7,6.3,6.7,7.2,6.2,6.1,6.4,7.2,7.4,7.9,6.4,6.3,6.1,7.7,6.3,6.4,6,6.9,6.7,6.9,5.8,6.8,6.7,6.7,6.3,6.5,6.2,5.9],"y":[6,5.1,5.9,5.6,5.8,6.6,4.5,6.3,5.8,6.1,5.1,5.3,5.5,5,5.1,5.3,5.5,6.7,6.9,5,5.7,4.9,6.7,4.9,5.7,6,4.8,4.9,5.6,5.8,6.1,6.4,5.6,5.1,5.6,6.1,5.6,5.5,4.8,5.4,5.6,5.1,5.1,5.9,5.7,5.2,5,5.2,5.4,5.1],"mode":"markers","name":"virginica","marker":{"color":"#66C2A5"}},{"type":"scatter","inherit":false,"x":[7,6.4,6.9,5.5,6.5,5.7,6.3,4.9,6.6,5.2,5,5.9,6,6.1,5.6,6.7,5.6,5.8,6.2,5.6,5.9,6.1,6.3,6.1,6.4,6.6,6.8,6.7,6,5.7,5.5,5.5,5.8,6,5.4,6,6.7,6.3,5.6,5.5,5.5,6.1,5.8,5,5.6,5.7,5.7,6.2,5.1,5.7],"y":[4.7,4.5,4.9,4,4.6,4.5,4.7,3.3,4.6,3.9,3.5,4.2,4,4.7,3.6,4.4,4.5,4.1,4.5,3.9,4.8,4,4.9,4.7,4.3,4.4,4.8,5,4.5,3.5,3.8,3.7,3.9,5.1,4.5,4.5,4.7,4.4,4.1,4,4.4,4.6,4,3.3,4.2,4.2,4.2,4.3,3,4.1],"mode":"markers","name":"versicolor","marker":{"color":"#FC8D62"}},{"type":"scatter","inherit":false,"x":[5.1,4.9,4.7,4.6,5,5.4,4.6,5,4.4,4.9,5.4,4.8,4.8,4.3,5.8,5.7,5.4,5.1,5.7,5.1,5.4,5.1,4.6,5.1,4.8,5,5,5.2,5.2,4.7,4.8,5.4,5.2,5.5,4.9,5,5.5,4.9,4.4,5.1,5,4.5,4.4,5,5.1,4.8,5.1,4.6,5.3,5],"y":[1.4,1.4,1.3,1.5,1.4,1.7,1.4,1.5,1.4,1.5,1.5,1.6,1.4,1.1,1.2,1.5,1.3,1.4,1.7,1.5,1.7,1.5,1,1.7,1.9,1.6,1.6,1.5,1.4,1.6,1.6,1.5,1.5,1.4,1.5,1.2,1.3,1.4,1.3,1.5,1.3,1.3,1.3,1.6,1.9,1.4,1.6,1.4,1.5,1.4],"mode":"markers","name":"setosa","marker":{"color":"#8DA0CB"}}],"layout":{"xaxis":{"title":"Sepal.Length"},"yaxis":{"title":"Petal.Length"},"hovermode":"closest","margin":{"b":40,"l":60,"t":25,"r":10}},"url":null,"width":null,"height":null,"source":"A","config":{"modeBarButtonsToRemove":["sendDataToCloud"]},"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->


# System Information


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
## [1] plotly_3.6.0     ggplot2_2.1.0    calibrate_1.7.2  MASS_7.3-45     
## [5] data.table_1.9.6
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.6        knitr_1.13         magrittr_1.5      
##  [4] munsell_0.4.3      colorspace_1.2-6   R6_2.1.2          
##  [7] httr_1.2.1         stringr_1.0.0      plyr_1.8.4        
## [10] tools_3.3.1        grid_3.3.1         gtable_0.2.0      
## [13] htmltools_0.3.5    assertthat_0.1     yaml_2.1.13       
## [16] digest_0.6.10      tibble_1.1         gridExtra_2.2.1   
## [19] RColorBrewer_1.1-2 tidyr_0.6.0        formatR_1.4       
## [22] viridis_0.3.4      base64enc_0.1-3    htmlwidgets_0.7   
## [25] evaluate_0.9       rmarkdown_1.0      stringi_1.1.1     
## [28] scales_0.4.0       jsonlite_1.0       chron_2.3-47
```

```r
# save.image(compress = TRUE, )
```
