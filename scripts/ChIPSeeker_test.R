#!/usr/bin/env Rscript

# this script will run all the sample workflow commands from the vignette, for testing:
# https://www.bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html

# Check for required packages
# and install
for (package in c("ChIPseeker","TxDb.Hsapiens.UCSC.hg19.knownGene","clusterProfiler","org.Hs.eg.db")) {
  # if(package %in% rownames(installed.packages()) == FALSE){install.packages(package, repos="http://cran.us.r-project.org")}
  if(package %in% rownames(installed.packages()) == FALSE){
    source("https://bioconductor.org/biocLite.R")
    biocLite(package)
    }
}

# source("https://bioconductor.org/biocLite.R")
# biocLite("clusterProfiler")

## loading packages
require(ChIPseeker)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
require(clusterProfiler)
library("org.Hs.eg.db")

plot_outdir <- "./"

# 4 ChIP profiling

# use built in sample files
files <- getSampleFiles()

# or using my own files:
# files <- c() # character vector of files paths
# convert to named list:
# files <- as.list(setNames(files,sapply(files,function(x) basename(dirname(x)))))

print(files)

peak <- readPeakFile(files[[4]])
peak

# 4.1 ChIP peaks coverage plot
pdf(paste0(plot_outdir, "/covplot.pdf"))
covplot(peak, weightCol="V5")
dev.off()

# 4.2 Profile of ChIP peaks binding to TSS regions

## promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
## tagMatrix <- getTagMatrix(peak, windows=promoter)
##
## to speed up the compilation of this vignettes, we use a precalculated tagMatrix
data("tagMatrixList")
tagMatrix <- tagMatrixList[[4]]

# 4.2.1 Heatmap of ChIP binding to TSS regions
pdf(paste0(plot_outdir, "/tagHeatmap.pdf"))
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
dev.off()

pdf(paste0(plot_outdir, "/peakHeatmap.pdf"))
peakHeatmap(files[[4]], TxDb=txdb, upstream=3000, downstream=3000, color="red")
dev.off()


# 4.2.2 Average Profile of ChIP peaks binding to TSS region
pdf(paste0(plot_outdir, "/AvgProfileTSSbingind.pdf"))
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), 
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.off()

pdf(paste0(plot_outdir, "/AvgProfileTSSbingind2.pdf"))
plotAvgProf2(files[[4]], TxDb=txdb, upstream=3000, downstream=3000, 
             xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.off()

pdf(paste0(plot_outdir, "/AvgProfileTSSbingind_conf.pdf"))
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 500)
dev.off()

# 4.3 Profile of ChIP peaks binding to start site of Exon/Intron

# 5 Peak Annotation
peakAnno <- annotatePeak(files[[4]], tssRegion=c(-3000, 3000), 
                         TxDb=txdb, annoDb="org.Hs.eg.db")

# 5.1 Visualize Genomic Annotation
pdf(paste0(plot_outdir, "/AnnoPie.pdf"))
plotAnnoPie(peakAnno)
dev.off()

pdf(paste0(plot_outdir, "/AnnoBar.pdf"))
plotAnnoBar(peakAnno)
dev.off()

pdf(paste0(plot_outdir, "/AnnoVennPie.pdf"))
vennpie(peakAnno)
dev.off()

pdf(paste0(plot_outdir, "/AnnoVennPie_UpSet.pdf"))
upsetplot(peakAnno, vennpie=TRUE)
dev.off()

# 5.2 Visualize distribution of TF-binding loci relative to TSS

pdf(paste0(plot_outdir, "/DistToTSS.pdf"))
plotDistToTSS(peakAnno, 
              title="Distribution of transcription factor-binding loci\nrelative to TSS")
dev.off()

# 6 Functional enrichment analysis

bp <- enrichGO(as.data.frame(peakAnno)$geneId, OrgDb='org.Hs.eg.db', ont="BP", readable=TRUE)
head(summary(bp), n=3)

# 7 ChIP peak data set comparison
# 7.1 Profile of several ChIP peak data binding to TSS region
# 7.1.1 Average profiles

## promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
## tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
## 
## to speed up the compilation of this vigenette, we load a precaculated tagMatrixList
data("tagMatrixList")

pdf(paste0(plot_outdir, "/AvgProf_default_tagMatrix.pdf"))
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
dev.off()


## resample = 500 by default, here use 100 to speed up the compilation of this vignette.
pdf(paste0(plot_outdir, "/AvgProf_default_tagMatrix_resample.pdf"))
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=100, facet="row")
dev.off()

# 7.1.2 Peak heatmaps
pdf(paste0(plot_outdir, "/tagHeatmap.pdf"))
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)
dev.off()

# 7.2 ChIP peak annotation comparision
# need a named list.. 
# as.list(setNames(rep(as.numeric(files[1]), length(files) - 1), files[-1]))
# files <- as.list(setNames(files,files))

# reset the names for the items
names(peakAnnoList) <- names(files)
peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-3000, 3000), verbose=FALSE)

pdf(paste0(plot_outdir, "/AnnoBar_list.pdf"))
plotAnnoBar(peakAnnoList,ylab=names(files))
dev.off()

pdf(paste0(plot_outdir, "/AnnoBar_dist_to_TSS.pdf"))
plotDistToTSS(peakAnnoList)
dev.off()


# 7.3 Functional profiles comparison
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster   = genes, 
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05, 
                           pAdjustMethod = "BH")
pdf(paste0(plot_outdir, "/compKEGG.pdf"))
plot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")
dev.off()


# 7.4 Overlap of peaks and annotated genes
genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)

pdf(paste0(plot_outdir, "/vennplot.pdf"))
vennplot(genes)
dev.off()


# 8 Statistical testing of ChIP seq overlap
# 8.1 Shuffle genome coordination
p <- GRanges(seqnames=c("chr1", "chr3"), 
             ranges=IRanges(start=c(1, 100), end=c(50, 130)))
shuffle(p, TxDb=txdb)


# 8.2 Peak overlap enrichment analysis
enrichPeakOverlap(queryPeak     = files[[4]], 
                  targetPeak    = unlist(files[1:4]), 
                  TxDb          = txdb, 
                  pAdjustMethod = "BH", 
                  nShuffle      = 50, 
                  chainFile     = NULL,
                  verbose       = FALSE)

# 9 Data Mining with ChIP seq data deposited in GEO
# 9.1 GEO data collection
getGEOspecies()
getGEOgenomeVersion()

hg19 <- getGEOInfo(genome="hg19", simplify=TRUE)
head(hg19)

# 9.2 Download GEO ChIP data sets
downloadGEObedFiles(genome="hg19", destDir="hg19")

gsm <- hg19$gsm[sample(nrow(hg19), 10)]
# downloadGSMbedFiles(gsm, destDir="hg19")


