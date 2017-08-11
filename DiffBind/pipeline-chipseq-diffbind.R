#!/usr/bin/env Rscript


# output width
options(width=300)

# load libraries
library(DiffBind, quietly=T)
library(RColorBrewer, quietly=T)
library(biomaRt, quietly=T)
library(ChIPpeakAnno, quietly=T)

# process command-line arguments (only arguments after --args)
args = commandArgs(trailingOnly=T)
outDir = args[1]
sampleSheetCsv = args[2]
genome = args[3] # hg19 or mm10
# blockFactor = args[4] 



# more functions
mycat <- function(text){cat(gsub(pattern = "\n", replacement = "  \n", x = text))}
sysinfo <- function(){
    # print custom information about the system
    
    # check if 'mycat' is loaded in case I copy/pasted this from elsewhere
    if( ! exists('mycat')) mycat <- function(text){cat(gsub(pattern = "\n", replacement = "  \n", x = text))}
    
    # system info for use on Linux with GNU tools installed
    # mycat(sprintf("System:\n%s\n%s", system("hostname", intern = TRUE), system("uname -srv", intern = TRUE)))
    # mycat(sprintf("System user:\n%s", system("whoami", intern = TRUE)))
    
    # dir
    # mycat(sprintf("System location:\n%s", system('pwd', intern = T, ignore.stderr = TRUE)))
    mycat(sprintf("System location:\n%s", getwd()))
    
    
    # repo info
    mycat(sprintf("Git Remote:\n%s\n", system('git remote -v', intern=T)))
    mycat(sprintf("Git branch and commit\n%s", system('printf "%s: %s" "$(git rev-parse --abbrev-ref HEAD)" "$(git rev-parse HEAD)"', 
                                                      intern = TRUE,  ignore.stderr = TRUE)))
    
    # date time
    mycat(sprintf("Time and Date of report creation:\n%s", system("date", intern = TRUE)))
    
    # R system info, packages, etc
    print(sessionInfo())
    print(Sys.info())
}

mytime <- function(){
    if( ! exists('mycat')) mycat <- function(text){cat(gsub(pattern = "\n", replacement = "  \n", x = text))}
    mycat(sprintf("Time: %s", date()))
}


mycat("Starting script")
mytime()
save.image(file = file.path(outDir, "loaded_args.Rdata"), compress = TRUE)




# extract differentially bound peaks and annotate them
generateDiffBindReport = function(dba, contrast, th=0.05, method=DBA_DESEQ2, reps=F, tss, mart.df, out.dir=".")
{
	library(DiffBind, quietly=T)
	library(ChIPpeakAnno, quietly=T)

	# extract contrast group names
	contrasts = as.matrix(dba.show(dba, bContrasts=T))
	group1 = contrasts[as.character(contrast), "Group1"]
	group2 = contrasts[as.character(contrast), "Group2"]

	# generate report (GRanges object)
	message("[generateDiffBindReport] generate diffbind report")
	message("[generateDiffBindReport] contrast num: ", contrast)
	message("[generateDiffBindReport] group1: ", group1)
	message("[generateDiffBindReport] group2: ", group2)
	message("[generateDiffBindReport] threshold: ", th)
	db.gr = dba.report(dba, contrast=contrast, method=method, th=th, bCounts=reps, DataType=DBA_DATA_GRANGES)

	message("[generateDiffBindReport] num sig peaks: ", length(db.gr))

	# annotate if any significant results were found
	if (length(db.gr)) {

		message("[generateDiffBindReport] annotate peaks")
		db.ann.gr = annotatePeakInBatch(db.gr, AnnotationData=tss, PeakLocForDistance="middle", FeatureLocForDistance="TSS", output="shortestDistance", multiple=T)

		# add gene symbols
		message("[generateDiffBindReport] add gene symbols")
		db.ann.df = merge(as.data.frame(db.ann.gr) , mart.df , by.x=c("feature"), by.y=c("ensembl_gene_id") , all.x=T)

		# keep just the relevant columns
		db.ann.df = db.ann.df[c(
			"seqnames","start","end","feature","external_gene_name","gene_biotype",
			"start_position","end_position",
			"insideFeature","distancetoFeature","shortestDistance","fromOverlappingOrNearest")]

		# merge bed and annotations
		ann.merged = merge(as.data.frame(db.gr), db.ann.df, by.x=c("seqnames","start","end"), by.y=c("seqnames","start","end") , all=T)

		# sort
		ann.merged$seqnames = as.character(ann.merged$seqnames)
		ann.merged = ann.merged[with(ann.merged, order(seqnames, start)), ]

		message("[generateDiffBindReport] save file")
		# generate file name
		th = format(th, nsmall=2)
		contrast.name = paste(group1, "-vs-", group2, sep="")
		if("Block1Val" %in% colnames(contrasts))
		{
			contrast.name = paste(contrast.name, ".blocking", sep="")
		}
		filename = paste(out.dir, "/diff_bind.", contrast.name, ".p", gsub(pattern="\\.", replacement="", x=th), ".csv", sep="")
		message("[generateDiffBindReport] save as: ", filename)
		write.csv(ann.merged, row.names=F, file=filename)
	}
}

message(" ========== load data ========== ")

# load data
sampleSheet = read.csv(sampleSheetCsv)
print(sampleSheet)
db = dba(sampleSheet=sampleSheet, bCorPlot=F, config=data.frame(AnalysisMethod=DBA_DESEQ2, th=0.05, cores=4))
print(db)
save.image(file = file.path(outDir, "loaded_args_and_data.Rdata"), compress = TRUE)


message(" ========== calculate binding matrix ========== ")

# calculate a binding matrix with scores based on read counts for every sample (affinity scores),
# rather than confidence scores for only those peaks called in a specific sample (occupancy scores)
db = dba.count(db, score=DBA_SCORE_RPKM)

# automatic contrasts
db = dba.contrast(db, categories=DBA_CONDITION, minMembers=2)

# save
dbContrastRData = paste(outDir, "/db.contrast.RData", sep="")
save(db, file=dbContrastRData)

message(" ========== generate plots ========== ")

# colors (combine multiple to prevent running out of colors)
colors = c(brewer.pal(9, "Set1"), brewer.pal(8, "Accent"), brewer.pal(8, "Dark2"))
colors = unique(colors)

# heatmap
heatmapPDF = paste(outDir, "/plot.heatmap.fpkm.pdf", sep="")
pdf(heatmapPDF, family="Palatino", pointsize=10)
dba.plotHeatmap(db, score=DBA_SCORE_RPKM, colScheme="Blues", colSideCols=colors, rowSideCols=colSideCols, cexRow=0.8, cexCol=0.8)
dev.off()

# pca based on conditions
pcaPDF = paste(outDir, "/plot.pca.fpkm.pdf", sep="")
pdf(pcaPDF, family="Palatino", pointsize=10)
dba.plotPCA(db, DBA_CONDITION, score=DBA_SCORE_RPKM, label=DBA_ID, vColors=colors, labelSize=0.6)
dev.off()

message(" ========== perform differential binding analysis ========== ")

# perform differential binding affinity analysis
db = dba.analyze(db, bFullLibrarySize=T, bCorPlot=F)

# save
dbAnalyzeRData = paste(outDir, "/db.analyze.RData", sep="")
save(db, file=dbAnalyzeRData)

# show contrasts
contrasts = dba.show(db, bContrasts=T)
contrasts = as.data.frame(contrasts)
print(contrasts)

message(" ========== retrieve biomart annotations ========== ")

# retrieve annotations
if (genome == "hg19") {
	martEns = useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", verbose=F)
}
if (genome == "mm10") {
	martEns = useMart(host="useast.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", verbose=F)
}
martEnsDF = getBM(attributes=c("ensembl_gene_id", "external_gene_name", "gene_biotype"), mart=martEns)
martEnsTSS = getAnnotation(mart=martEns, featureType="TSS")

message(" ========== generate reports ========== ")

# generate report for each possible contrast with different cutoffs
for (i in 1:length(row.names(contrasts))) {
	print(contrasts[i,])
	# using try to prevent "execution halted" that kills script if there are no significant results
	try(generateDiffBindReport(dba=db, contrast=i, th=1.00, method=DBA_DESEQ2, tss=martEnsTSS, mart.df=martEnsDF, reps=T, out.dir=outDir))
	try(generateDiffBindReport(dba=db, contrast=i, th=0.20, method=DBA_DESEQ2, tss=martEnsTSS, mart.df=martEnsDF, reps=T, out.dir=outDir))
	try(generateDiffBindReport(dba=db, contrast=i, th=0.05, method=DBA_DESEQ2, tss=martEnsTSS, mart.df=martEnsDF, reps=T, out.dir=outDir))
}

# repeat with blocking factor if blocking factor parameter was passed
if (!is.na(args[4])) {
	message(" ========== re-analyze with blocking factor ========== ")

	# automatic contrasts with blocking factor
	db = dba.contrast(db, categories=DBA_CONDITION, block=DBA_REPLICATE)
	db = dba.analyze(db, bFullLibrarySize=T, bCorPlot=F)
	dbAnalyzeBlockingRData = paste(outDir, "/db.analyze.blocking.RData", sep="")
	save(db, file=dbAnalyzeBlockingRData)

	# show contrasts
	contrasts = dba.show(db, bContrasts=T)
	contrasts = as.data.frame(contrasts)
	print(contrasts)

	message(" ========== generate reports with blocking factor ========== ")

	# generate report for each possible contrast
	for (i in 1:length(row.names(contrasts))) {
		print(contrasts[i,])
		# using try to prevent "execution halted" that kills script if there are no significant results
		try(generateDiffBindReport(dba=db, contrast=i, th=1.00, method=DBA_DESEQ2_BLOCK, tss=martEnsTSS, mart.df=martEnsDF, reps=T, out.dir=outDir))
		try(generateDiffBindReport(dba=db, contrast=i, th=0.20, method=DBA_DESEQ2_BLOCK, tss=martEnsTSS, mart.df=martEnsDF, reps=T, out.dir=outDir))
		try(generateDiffBindReport(dba=db, contrast=i, th=0.05, method=DBA_DESEQ2_BLOCK, tss=martEnsTSS, mart.df=martEnsDF, reps=T, out.dir=outDir))
	}
}




mycat("Finished script")
mytime()
sysinfo()
save.image(file = file.path(outDir, "final.Rdata"), compress = TRUE)
# end
