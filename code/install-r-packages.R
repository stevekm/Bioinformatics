#!/usr/bin/env Rscript

install_bioconductor <- "FALSE"
install_baseR <- "FALSE"

# get commands passed to the script
args <- commandArgs(TRUE)

# test if there is at least one argument
if (length(args)==0) {
    args[1] <- ""
} 

if(as.character(args[1]) == "base"){
    install_baseR <- "TRUE"
} else if(args[1] == "bioc"){
    install_bioconductor <- "TRUE"
} else if(args[1] == "all"){
    install_baseR <- "TRUE"
    install_bioconductor <- "TRUE"
} else {
    cat("No pacakges were installed. Please specify one of the following:\nbase\nbioc\nall\n\n")
}

if(install_bioconductor){
    write("Checking installation of Bioconductor packages", stderr())
    # new version:
    # Bioconductor packages
    for (package in c("ChIPseeker","TxDb.Hsapiens.UCSC.hg19.knownGene","clusterProfiler","org.Hs.eg.db","wavClusteR","DiffBind",
                      "biomaRt","ChIPpeakAnno", 'preprocessCore')) {
        # if(package %in% rownames(installed.packages()) == FALSE){install.packages(package, repos="http://cran.us.r-project.org")}
        if(package %in% rownames(installed.packages()) == FALSE){
            source("http://bioconductor.org/biocLite.R") # https:// URLs are not supported
            biocLite(package,suppressUpdates=FALSE,ask = FALSE) # don't update dependency packages, dont ask
        }
    }
}


if(install_baseR){
    write("Checking installation of base R packages", stderr())
    package_list<-c("ggplot2", "grid", "plyr", "knitr", "VennDiagram", 
                    "gridExtra", "datasets", "digest", "Hmisc", "xtable", "reshape2", 
                    "data.table", "scales", "corrplot", "RColorBrewer", "lattice", 
                    "gplots", "MASS", "stringr", "flsa", "genlasso", "optparse", 
                    "pastecs", "plotrix", "zoo", "reshape", "chron","UpSetR", "plotly") 
    
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
                     warning = function(e){write(paste0("Package not installed: ", p), stderr())},
                     error = function(e){write(paste0("Package not installed: ", p), stderr())})
            
            # try to install & load the packages, skip to next loop iteration upon failure
            tryCatch(install.packages(p,repos="http://cran.rstudio.com/"),warning = next)
            tryCatch(library(p,character.only=TRUE,verbose=FALSE),warning = next)
        }
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
