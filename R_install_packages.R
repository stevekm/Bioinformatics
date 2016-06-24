#!/usr/bin/env Rscript

# this script will install Bioconductor and R packages

# Bioconductor packages
for (package in c("ChIPseeker","TxDb.Hsapiens.UCSC.hg19.knownGene","clusterProfiler","org.Hs.eg.db","wavClusteR","DiffBind",
                  "biomaRt","ChIPpeakAnno")) {
  # if(package %in% rownames(installed.packages()) == FALSE){install.packages(package, repos="http://cran.us.r-project.org")}
  if(package %in% rownames(installed.packages()) == FALSE){
    source("https://bioconductor.org/biocLite.R")
    biocLite(package)
  }
}

# R packages
package_list <-c("UpSetR","ggplot2","grid","plyr","knitr","VennDiagram","gridExtra","datasets","digest","Hmisc","xtable",
                 "reshape2","data.table","scales","corrplot","RColorBrewer","lattice","gplots","MASS","stringr",
                 "flsa","genlasso","ggplot2","optparse","pastecs","plotrix","reshape2","zoo","optparse","gplots","ggplot2","reshape",
                 'optparse','MASS') # 'preprocessCore'
for (p in package_list) 
  if (!require(p,character.only=TRUE,quietly=TRUE,warn.conflicts=FALSE)) {
    install.packages(p,repos="http://cran.rstudio.com/") 
    library(p,character.only=TRUE,verbose=FALSE)
  }

