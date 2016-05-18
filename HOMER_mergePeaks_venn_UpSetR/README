This simple workflow uses the venn.txt output from HOMER's mergePeaks command to visualize the overlaps between different sets of peaks (such as .bed files from a ChIP-Seq experiment), using the UpSetR package for R version 3.3.0. 

In this example, the `venn.txt` file would have been created by using a HOMER command such as this:

```

module load homer/v4.6
mergePeaks H3K27AC.bed H3K27ME3.bed H3K4ME3.bed H3K9AC.bed gencode.bed -prefix mergepeaks -venn venn.txt -matrix matrix.txt


```

Pass the `venn.txt` file to the `multi_peaks_UpSet_plot.R` script like this:

```

Rscript --vanilla multi_peaks_UpSet_plot.R "SampleID" venn.txt

```
