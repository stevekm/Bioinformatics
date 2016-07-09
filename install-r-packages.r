#!/usr/bin/env Rscript

write("Checking installation of R packages", stderr())

# new version:
# Bioconductor packages
for (package in c("ChIPseeker","TxDb.Hsapiens.UCSC.hg19.knownGene","clusterProfiler","org.Hs.eg.db","wavClusteR","DiffBind",
                  "biomaRt","ChIPpeakAnno")) {
  # if(package %in% rownames(installed.packages()) == FALSE){install.packages(package, repos="http://cran.us.r-project.org")}
  if(package %in% rownames(installed.packages()) == FALSE){
    source("http://bioconductor.org/biocLite.R") # https:// URLs are not supported
    biocLite(package,suppressUpdates=FALSE,ask = FALSE) # don't update dependency packages, dont ask
  }
}

package_list<-c("ggplot2", "grid", "plyr", "knitr", "VennDiagram", 
                "gridExtra", "datasets", "digest", "Hmisc", "xtable", "reshape2", 
                "data.table", "scales", "corrplot", "RColorBrewer", "lattice", 
                "gplots", "MASS", "stringr", "flsa", "genlasso", "optparse", 
                "pastecs", "plotrix", "zoo", "reshape", "chron","UpSetR",'preprocessCore') 

for(p in package_list){
  # check if the package is already installed
  if(require(p,character.only=TRUE,quietly=TRUE,warn.conflicts=FALSE)){
    write(paste0("Package already installed: ", p), stderr())
  }
  
  # check if package can't be loaded
  if(!require(p,character.only=TRUE,quietly=TRUE,warn.conflicts=FALSE)){
    write(paste0("Attempting to install package: ",p), stderr())
    # try to install & load the packages, give a message upon failure
    tryCatch(install.packages(p,repos="http://cran.rstudio.com/"),
             warning = function(e){write(paste0("Failed to install pacakge: ", p), stderr())},
             error = function(e){write(paste0("Failed to install pacakge: ", p), stderr())})
    tryCatch(library(p,character.only=TRUE,verbose=FALSE),
             warning = function(e){write(paste0("Failed to install pacakge: ", p), stderr())},
             error = function(e){write(paste0("Failed to install pacakge: ", p), stderr())})
    
    # try to install & load the packages, skip to next loop iteration upon failure
    tryCatch(install.packages(p,repos="http://cran.rstudio.com/"),warning = next)
    tryCatch(library(p,character.only=TRUE,verbose=FALSE),warning = next)
  }
}


# old version :
# packages = c("flsa", "genlasso", "optparse", "ggplot2", "pastecs", "plotrix",
#              "reshape2", "zoo", "plyr", "data.table", "gridExtra", "scales", "grid",
#              "RColorBrewer", "stringr", "chron", "corrplot")
# 
# for (p in packages)
#   if (!suppressPackageStartupMessages(require(p, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE))) {
#     install.packages(p, repos="http://cran.us.r-project.org")
#   }
