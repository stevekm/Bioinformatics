#!/usr/bin/env Rscript

## USAGE: ....R /path/to/peaks.bed /path/to/outdir my_sample_ID
## DESCRIPTION:

# get script args
args <- commandArgs(TRUE)

cat("\nScript args are:\n")
print(args)
input_peaks_file <- args[1]
output_directory <- args[2]
sampleID <- args[3]

# promoter_proximal = 3000 # Extending promoter upstream and downstream by nt
# Rscript --vanilla code/chipseq-peakanno.r -g $genome -d $promoter_proximal -o $outdir $peaks
cat("\nLoading packages...\n")
library("ChIPseeker")
library("clusterProfiler")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")

cat("\nLoading hg19 database...\n")
txdb <- get("TxDb.Hsapiens.UCSC.hg19.knownGene")
promoter_dist <- 3000

cat("\nReading peaks file...\n")
peak <- readPeakFile(input_peaks_file)

cat("\nMaking Chrom Coverages plot...\n")
peaks_coverage_plot_file <- file.path(output_directory, "peaks-coverage.pdf")
sample_title <- paste0(sampleID, " ChIP Peaks over Chromosomes")
pdf(file = peaks_coverage_plot_file)
covplot(peak, weightCol="V5", title = sample_title) # title = "ChIP Peaks over Chromosomes"
dev.off()

cat("\nGetting peak annotations...\n")
peakAnno <- annotatePeak(peak, tssRegion=c(-promoter_dist, promoter_dist), 
                         TxDb=txdb, 
                         annoDb="org.Hs.eg.db")

cat("\nMaking Peak Anno pie chart...\n")
anno_piechart_plot_file <- file.path(output_directory, "anno-piechart.pdf")
sample_title <- paste0("\n\n", sampleID, " Peak Types")
pdf(file = anno_piechart_plot_file, height = 8, width = 8)
plotAnnoPie(peakAnno, main = sample_title)
dev.off()

cat("\nMaking Upset plot...\n")
upset_plot_file <- file.path(output_directory, "upsetplot.pdf")
sample_title <- paste0(sampleID, " Peak Overlaps")
pdf(file = upset_plot_file, width = 9, height = 4.5, onefile = F)
upsetplot(peakAnno,vennpie=TRUE) 
text(x = 0, y = 1, sample_title) # add a title
dev.off()

cat("\nSaving table...\n")
peak_anno_table_file <- file.path(output_directory, "peak_anno.tsv")
write.table(peakAnno, quote=FALSE, sep="\t", row.names =FALSE, file=peak_anno_table_file)

sessionInfo()
