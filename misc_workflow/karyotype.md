[[source](https://www.biostars.org/p/269857/)]
# Creating chromosome karyotype plot with R and ggplot2

There are numerous resources for creating karyotype and ideogram plots, such as those posted [here](https://www.biostars.org/p/378/). You can use a variety of software packages for this, along with various R Bioconductor packages. 

But what happens when you want to customize your plots beyond what is made easily available by your library or software package? Or you are stuck in 'dependency hell' and have issues integrating those software packages and libraries into your automated workflow? You might consider just drawing the plots yourself instead. 

In this example, I will show the steps you can take to create some basic versions of chromosome plots, which you can easily customize to your liking due to the fact that it is drawn completely in R using `ggplot2`. I will also be using the `ggrepel` package, which helps with placing non-overlapping text labels on the plot, and the `scales` package which helps to format values displayed on the plot axis.

# Load Data

I have included the R code here necessary to re-create the sample datasets directly as they were read in from files. You can replace this step with the code used to read your data from file (e.g. `read.delim()`). In this example I am using some copy number alteration data generated with CNVkit, and chromosome & centromere values taken from [here](https://www.biostars.org/p/2349/)

```r
library("ggplot2") # for the plot
library("ggrepel") # for spreading text labels on the plot
library("scales") # for axis labels notation

# insert your steps to load data from tabular files or other sources here; 
# dummy datasets taken directly from files shown in this example

# data with the copy number alterations for the sample
sample_cns <- structure(list(gene = c("AC116366.7", "ANKRD24", "APC", "SNAPC3", 
"ARID1A", "ATM", "BOD1L1", "BRCA1", "C11orf65", "CHD5", "COL1A1", 
"CTC-554D6.1", "CTD-2047H16.4", "DNAH9", "DOT1L", "DPYD", "DPYD-AS1", 
"EP300", "EP300-AS1", "ERBB2", "FGFR2", "KMT2A", "KMT2B", "LRP1B", 
"MAP3K1"), chromosome = c("chr5", "chr19", "chr5", "chr9", "chr1", 
"chr11", "chr4", "chr17", "chr11", "chr1", "chr17", "chr5", "chr17", 
"chr17", "chr19", "chr1", "chr1", "chr22", "chr22", "chr17", 
"chr10", "chr11", "chr19", "chr2", "chr5"), start = c(131893016L, 
4183350L, 112043414L, 15465517L, 27022894L, 108098351L, 13571634L, 
41197694L, 108180886L, 6166339L, 48262862L, 112162804L, 78326739L, 
11501815L, 2164183L, 97544531L, 97564044L, 41489008L, 41572250L, 
37844086L, 123239094L, 118307227L, 36208920L, 140990754L, 56111400L
), end = c(131978056L, 4224502L, 112179823L, 15465578L, 27107247L, 
108236235L, 13629211L, 41276113L, 108236235L, 6240083L, 48278874L, 
112179823L, 78367298L, 11872844L, 2229791L, 98386478L, 97771853L, 
41574960L, 41574960L, 37884297L, 123353331L, 118392887L, 36229458L, 
142888298L, 56189507L), log2 = c(-0.850333, -0.802459, -0.850333, 
1.68765, -0.828046, -0.883559, 0.495105, 0.51503, -0.883559, 
-0.828046, 0.51503, -0.850333, 0.51503, -0.801607, -0.802459, 
-0.828046, -0.828046, 0.517179, 0.517179, 0.51503, -0.865372, 
-0.883559, 0.970523, -0.809056, -0.850333), depth = c(823.473, 
240.685, 723.721, 4325.57, 596.063, 560.472, 1563.44, 1609.96, 
703.526, 411.25, 1586.75, 986.95, 2779.47, 643.981, 219.58, 654.69, 
648.597, 1488.49, 2631.62, 2144.13, 893.806, 222.718, 985.121, 
1112.21, 571.51), weight = c(16.7856, 17.0764, 31.7769, 0.557449, 
23.296, 39.5052, 34.3571, 23.5551, 15.7455, 25.5399, 28.9927, 
22.9053, 23.2428, 52.522, 26.4509, 18.9309, 3.71943, 27.8139, 
7.18582, 18.225, 21.8383, 43.5557, 31.4704, 58.351, 19.5343), 
    cn = c(1L, 1L, 1L, 7L, 1L, 1L, 3L, 3L, 1L, 1L, 3L, 1L, 3L, 
    1L, 1L, 1L, 1L, 3L, 3L, 3L, 1L, 1L, 4L, 1L, 1L), probes = c(897L, 
    508L, 897L, 51L, 1052L, 434L, 370L, 847L, 434L, 1052L, 847L, 
    897L, 847L, 284L, 508L, 1052L, 1052L, 125L, 125L, 847L, 157L, 
    434L, 66L, 226L, 897L)), .Names = c("gene", "chromosome", 
"start", "end", "log2", "depth", "weight", "cn", "probes"), row.names = c(NA, 
25L), class = "data.frame")

# > head(sample_cns)
#         gene chromosome     start       end      log2    depth    weight cn probes
# 1 AC116366.7       chr5 131893016 131978056 -0.850333  823.473 16.785600  1    897
# 2    ANKRD24      chr19   4183350   4224502 -0.802459  240.685 17.076400  1    508
# 3        APC       chr5 112043414 112179823 -0.850333  723.721 31.776900  1    897
# 4     SNAPC3       chr9  15465517  15465578  1.687650 4325.570  0.557449  7     51
# 5     ARID1A       chr1  27022894  27107247 -0.828046  596.063 23.296000  1   1052
# 6        ATM      chr11 108098351 108236235 -0.883559  560.472 39.505200  1    434

# hg19 chromosome sizes
chrom_sizes <- structure(list(V1 = c("chrM", "chr1", "chr2", "chr3", "chr4", 
"chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
"chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
"chr20", "chr21", "chr22", "chrX", "chrY"), V2 = c(16571L, 249250621L, 
243199373L, 198022430L, 191154276L, 180915260L, 171115067L, 159138663L, 
146364022L, 141213431L, 135534747L, 135006516L, 133851895L, 115169878L, 
107349540L, 102531392L, 90354753L, 81195210L, 78077248L, 59128983L, 
63025520L, 48129895L, 51304566L, 155270560L, 59373566L)), .Names = c("V1", 
"V2"), class = "data.frame", row.names = c(NA, -25L))

# > head(chrom_sizes)
#     V1        V2
# 1 chrM     16571
# 2 chr1 249250621
# 3 chr2 243199373
# 4 chr3 198022430
# 5 chr4 191154276
# 6 chr5 180915260


# hg19 centromere locations
centromeres <- structure(list(X.bin = c(23L, 20L, 2L, 1L, 14L, 16L, 1L, 14L, 
1L, 1L, 10L, 1L, 15L, 13L, 1L, 1L, 11L, 13L, 1L, 1L, 1L, 12L, 
10L, 10L), chrom = c("chr1", "chr2", "chr3", "chr4", "chr5", 
"chr6", "chr7", "chr8", "chr9", "chrX", "chrY", "chr10", "chr11", 
"chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
"chr19", "chr20", "chr21", "chr22"), chromStart = c(121535434L, 
92326171L, 90504854L, 49660117L, 46405641L, 58830166L, 58054331L, 
43838887L, 47367679L, 58632012L, 10104553L, 39254935L, 51644205L, 
34856694L, 16000000L, 16000000L, 17000000L, 35335801L, 22263006L, 
15460898L, 24681782L, 26369569L, 11288129L, 13000000L), chromEnd = c(124535434L, 
95326171L, 93504854L, 52660117L, 49405641L, 61830166L, 61054331L, 
46838887L, 50367679L, 61632012L, 13104553L, 42254935L, 54644205L, 
37856694L, 19000000L, 19000000L, 20000000L, 38335801L, 25263006L, 
18460898L, 27681782L, 29369569L, 14288129L, 16000000L), ix = c(1270L, 
770L, 784L, 447L, 452L, 628L, 564L, 376L, 411L, 583L, 105L, 341L, 
447L, 304L, 3L, 3L, 3L, 354L, 192L, 125L, 410L, 275L, 22L, 3L
), n = c("N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", 
"N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N"
), size = c(3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 
3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 
3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 
3000000L, 3000000L, 3000000L, 3000000L, 3000000L), type = c("centromere", 
"centromere", "centromere", "centromere", "centromere", "centromere", 
"centromere", "centromere", "centromere", "centromere", "centromere", 
"centromere", "centromere", "centromere", "centromere", "centromere", 
"centromere", "centromere", "centromere", "centromere", "centromere", 
"centromere", "centromere", "centromere"), bridge = c("no", "no", 
"no", "no", "no", "no", "no", "no", "no", "no", "no", "no", "no", 
"no", "no", "no", "no", "no", "no", "no", "no", "no", "no", "no"
)), .Names = c("X.bin", "chrom", "chromStart", "chromEnd", "ix", 
"n", "size", "type", "bridge"), class = "data.frame", row.names = c(NA, 
-24L))

# > head(centromeres)
#   X.bin chrom chromStart  chromEnd   ix n    size       type bridge
# 1    23  chr1  121535434 124535434 1270 N 3000000 centromere     no
# 2    20  chr2   92326171  95326171  770 N 3000000 centromere     no
# 3     2  chr3   90504854  93504854  784 N 3000000 centromere     no
# 4     1  chr4   49660117  52660117  447 N 3000000 centromere     no
# 5    14  chr5   46405641  49405641  452 N 3000000 centromere     no
# 6    16  chr6   58830166  61830166  628 N 3000000 centromere     no

```

# Adjust Data

We need to take some steps in order to prepare our data for plotting. That includes making sure that column labels are consistent across datasets, and setting up an ordered factor level for the chromosomes. 

```r
# set the column names for the datasets
# IMPORTANT: fields common across datasets should have the same name in each
colnames(chrom_sizes) <- c("chromosome", "size")
colnames(centromeres) <- c('bin', "chromosome", 'start', 'end',
                       'ix', 'n', 'size', 'type', 'bridge')

# create an ordered factor level to use for the chromosomes in all the datasets
chrom_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                 "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
                 "chr22", "chrX", "chrY", "chrM")
chrom_key <- setNames(object = as.character(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
                                              12, 13, 14, 15, 16, 17, 18, 19, 20, 
                                              21, 22, 23, 24, 25)), 
                      nm = chrom_order)
chrom_order <- factor(x = chrom_order, levels = rev(chrom_order))

# convert the chromosome column in each dataset to the ordered factor
chrom_sizes[["chromosome"]] <- factor(x = chrom_sizes[["chromosome"]], 
                                      levels = chrom_order)
sample_cns[["chromosome"]] <- factor(x = sample_cns[["chromosome"]], 
                                     levels = chrom_order)
centromeres[["chromosome"]] <- factor(x = centromeres[["chromosome"]], 
                                      levels = chrom_order)
# mark 'gain' or 'loss' for each Copy Number Alteration (CNA) in the sample dataset
sample_cns$CNA <- ""
sample_cns$CNA[sample_cns$cn <=1] <- "loss"
sample_cns$CNA[sample_cns$cn >= 3] <- "gain"
# create a color key for the plot
group.colors <- c(gain = "red", loss = "blue")

```

# Make Plot

Instead of creating a bar plot off of the chromosome factor levels, we will instead convert the chromosome factors to their numeric value to allow for more fine-tuned spacing, and then display the original chromosome ID on the axis to give the appearance of a discrete scale. We will also use `geom_rect` to create the bars and bands in the plot, instead of something like `geom_bar`, so that we can have more control over the shapes being drawn. If you don't want to use `geom_text_repel`, you can substitute with `geom_text`. 

```r
ggplot(data = chrom_sizes) + 
    # base rectangles for the chroms, with numeric value for each chrom on the x-axis
    geom_rect(aes(xmin = as.numeric(chromosome) - 0.2, 
                  xmax = as.numeric(chromosome) + 0.2, 
                  ymax = size, ymin = 0), 
              colour="black", fill = "white") + 
    # rotate the plot 90 degrees
    coord_flip() +
    # black & white color theme 
    theme(axis.text.x = element_text(colour = "black"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank()) + 
    # give the appearance of a discrete axis with chrom labels
    scale_x_discrete(name = "chromosome", limits = names(chrom_key)) +
    # add bands for centromeres
    geom_rect(data = centromeres, aes(xmin = as.numeric(chromosome) - 0.2, 
                                      xmax = as.numeric(chromosome) + 0.2, 
                                      ymax = end, ymin = start)) +
    # add bands for CNA value
    geom_rect(data = sample_cns, aes(xmin = as.numeric(chromosome) - 0.2, 
                                     xmax = as.numeric(chromosome) + 0.2, 
                                     ymax = end, ymin = start, fill = CNA)) + 
    scale_fill_manual(values = group.colors) +
    # add 'gain' gene markers
    geom_text_repel(data = subset(sample_cns, sample_cns$CNA == "gain"), 
                    aes(x = chromosome, y = start, label = gene), 
                    color = "red", show.legend = FALSE) +
    # add 'loss' gene markers
    geom_text_repel(data = subset(sample_cns, sample_cns$CNA == "loss"), 
                    aes(x = chromosome, y = start, label = gene ), 
                    color = "blue", show.legend = FALSE) +
    ggtitle("Copy Number Alterations") +
    # supress scientific notation on the y-axis
    scale_y_continuous(labels = comma) +
    ylab("region (bp)")
```

# Results
![sample_chrom_plot](https://user-images.githubusercontent.com/10505524/29837704-e53fcb00-8cc6-11e7-9fe1-1d009abb2e0d.png)

  [1]: http://![sample_chrom_plot](https://user-images.githubusercontent.com/10505524/29837704-e53fcb00-8cc6-11e7-9fe1-1d009abb2e0d.png)

# Conclusion

This tutorial should help you create some basic chromosome plots, which you can further customize to suit your needs. 

# Extra Resources & Ideas

- [How to plot positions along a chromosome graphic](https://stackoverflow.com/questions/33727432/how-to-plot-positions-along-a-chromosome-graphic)

- [`ggplot2` reference manual](http://ggplot2.tidyverse.org/reference/)

- [`ggplot2` cheatsheet](http://www.rstudio.com/wp-content/uploads/2015/12/ggplot2-cheatsheet.pdf)

# Software

- R version 3.2.3

- `ggplot2_2.2.1`

- `ggrepel_0.6.5`

- `scales_0.4.1`
