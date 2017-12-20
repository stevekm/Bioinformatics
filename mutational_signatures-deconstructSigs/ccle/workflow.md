# CCLE mutation signatures

Running the [`deconstructSigs`](https://github.com/raerose01/deconstructSigs) package to find genomic mutational signatures in a publically available RNASeq tumor dataset from the Broad Institute here:

https://portals.broadinstitute.org/ccle/data

- using the file: `ccle2maf_081117.txt`

# Workflow

The workflow is divided into discrete steps to aid execution on the NYULMC phoenix compute cluster

```
# load and clean the CCLE data
$ Rscript ccle_clean.R

# generate signatures for all 1400 samples using the `ccle_make_signatures.R` script; 
# run on the compute cluster
$ ./qsub_ccle_make_signatures.sh

# extract the signature weights from all the samples
$ Rscript ccle_compare_signatures.R

```


# Software

Test with CentOS 6, R 3.2.3 / 3.3.0, and the following:

```
> sessionInfo()
R version 3.2.3 (2015-12-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS release 6.9 (Final)

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] BSgenome.Hsapiens.UCSC.hg19_1.4.0 BSgenome_1.38.0                   rtracklayer_1.30.4               
 [4] Biostrings_2.38.4                 XVector_0.10.0                    GenomicRanges_1.22.4             
 [7] GenomeInfoDb_1.6.3                IRanges_2.4.8                     S4Vectors_0.8.11                 
[10] BiocGenerics_0.16.1               deconstructSigs_1.8.0             reshape2_1.4.2                   
[13] ggplot2_2.2.1                    

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.12               magrittr_1.5               GenomicAlignments_1.6.3    zlibbioc_1.16.0           
 [5] BiocParallel_1.4.3         munsell_0.4.3              colorspace_1.3-2           rlang_0.1.1               
 [9] stringr_1.2.0              plyr_1.8.4                 tools_3.2.3                SummarizedExperiment_1.0.2
[13] grid_3.2.3                 Biobase_2.30.0             gtable_0.2.0               lambda.r_1.1.9            
[17] futile.logger_1.4.3        lazyeval_0.2.0             tibble_1.3.3               futile.options_1.0.0      
[21] bitops_1.0-6               RCurl_1.95-4.8             stringi_1.1.5              Rsamtools_1.22.0          
[25] scales_0.4.1               XML_3.98-1.9              
```
