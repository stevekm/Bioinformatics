#!/usr/bin/env Rscript

# write out the citations for a list of R packages
# this will give whatever is listed as the package citation, which will hopefully
# include a BibTeX format citation as well

# list of package names
package_list<-c("ggplot2", "grid", "plyr", "knitr", "VennDiagram", 
                "gridExtra", "datasets", "digest", "Hmisc", "xtable", "reshape2", 
                "data.table", "scales", "corrplot", "RColorBrewer", "lattice", 
                "gplots", "MASS", "stringr", "flsa", "genlasso", "optparse", 
                "pastecs", "plotrix", "zoo", "reshape", "chron","UpSetR",'preprocessCore',"ChIPseeker","TxDb.Hsapiens.UCSC.hg19.knownGene","clusterProfiler","org.Hs.eg.db","wavClusteR","DiffBind",
                "biomaRt","ChIPpeakAnno") 

# write them to the file...
sink("/ifs/home/kellys04/references.bib")
for(pkg in package_list){
    cat(paste0("-----------------\nR package:\n",pkg,"\n"))
    
    # only if they don't return an error...
    try(print(citation(pkg)),silent = TRUE)
    cat(paste0("-----------------\n"))
}
sink()
