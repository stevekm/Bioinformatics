#!/usr/bin/Rscript
# called by scripts-alignment-summary.sh
# USAGE: alignment_summary_stats.R /path/to/stats.tsv /path/to/outputdir
# create dual barplots to visualize alignment reads


# get commands passed to the script
args <- commandArgs()
# /path/to/stats.tsv will be args[6]
# /path/to/outputdir will be args[7] 

##################
# Set up the data

# get the reads.tsv from the script call
ReadsTable<-read.csv(args[6],row.names=1)

# calculate the percent alignment
ReadsTable$PercentAligned=signif(c(ReadsTable$Aligned.Reads / ReadsTable$Total.Reads)*100,digits = 4 )
# calculate the percent deduplicated
ReadsTable$PercentDedup=signif(c(ReadsTable$De.duplicated.alignments / ReadsTable$Total.Reads)*100,digits = 4 )
# calculate the number of unaligned reads
ReadsTable$UnalignedReads=c(ReadsTable$Total.Reads - ReadsTable$Aligned.Reads)
#calculate the percent unaligned reads
ReadsTable$PcntUnalignedReads=signif(c(ReadsTable$UnalignedReads / ReadsTable$Total.Reads)*100,digits=4)
# calculate the number of duplicated reads
ReadsTable$Duplicated=c(ReadsTable$Aligned.Reads - ReadsTable$De.duplicated.alignments)
# calculate the percent of duplicated
ReadsTable$PercentDup=signif(c(ReadsTable$Duplicated / ReadsTable$Total.Reads)*100,digits = 4 )

# save the values to be plotted into a transposed matrix, since thats what the barplot() likes
# # first just get the columns we want
Dup_Raw_Reads_df<-ReadsTable[,which(colnames(ReadsTable) %in% c("De.duplicated.alignments","Duplicated","UnalignedReads")) ] 
# reorder the columns because R is dumb
Dup_Raw_Reads_df<-Dup_Raw_Reads_df[c("De.duplicated.alignments","Duplicated","UnalignedReads")]

# Dup_Raw_Reads_Matrix<-t(as.matrix( ReadsTable[,which(colnames(ReadsTable) %in% c("De.duplicated.alignments","Duplicated","UnalignedReads")) ] ))
Dup_Raw_Reads_Matrix<-t(as.matrix(Dup_Raw_Reads_df))
# # divid the number of reads by 1million
Dup_Raw_Reads_Matrix<-signif(Dup_Raw_Reads_Matrix/1000000,digits = 4)

# # first just get the columns we want
Dup_Pcnt_Reads_df<-ReadsTable[,which(colnames(ReadsTable) %in% c("PercentDedup","PercentDup","PcntUnalignedReads")) ]
# reorder
Dup_Pcnt_Reads_df<-Dup_Pcnt_Reads_df[c("PercentDedup","PercentDup","PcntUnalignedReads")]
Dup_Pcnt_Reads_Matrix<-t(as.matrix(Dup_Pcnt_Reads_df))
# Dup_Pcnt_Reads_Matrix<-t(as.matrix( ReadsTable[,which(colnames(ReadsTable) %in% c("PercentDedup","PercentDup","PcntUnalignedReads")) ] ))

# Dup_Raw_Reads_Matrix
# Dup_Pcnt_Reads_Matrix
# 
# #####################
# # Set up the plots
# 
BARPLOT_COLORS<-c("blue","purple","red")

# setup the matrix for the plot layout
Raw_Reads_Matrix_matrix<-structure(c(1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 1L, 2L, 
                                     2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 
                                     3L, 3L, 3L, 3L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L),
                                   .Dim = c(11L,4L), 
                                   .Dimnames = list(NULL, c("V1", "V2", "V3", "V4")))

# setup the panel layout
# layout(Raw_Reads_Matrix_matrix)
# # layout.show(max(Raw_Reads_Matrix_matrix))  # use this to preview the layout
# 
# par(mar=c(0,0,4,0)) # need to set this for some reason
# # on mar: A numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot. 
# # The default is c(5, 4, 4, 2) + 0.1.
# 
# # call blank plot to fill the first panel
# plot(1,type='n',axes=FALSE,xlab="",ylab="",main = "Sequencing Reads",cex.main=2) 
# 
# # set up the Legend in the first panel
# legend("bottom",legend=c("Deduplicated","Duplicated","Unaligned"),fill=BARPLOT_COLORS,bty = "n",ncol=length(BARPLOT_COLORS),cex=1.0)
# 
# # set the mar for the remaining panels
# par(mar=c(6,5,0,4))
# # create barplot for the two matrices
# barplot(Dup_Raw_Reads_Matrix,horiz = T,col=BARPLOT_COLORS,border=NA,las=1,cex.names=1,xlab="Number of reads (millions)")
# barplot(Dup_Pcnt_Reads_Matrix,horiz = T,col=BARPLOT_COLORS,border=NA,las=1,cex.names=1,xlab="Percent of reads")

# # make a filename and filepath for the output
# OutFile_name<-paste0(dirname(args[7]),"_read_barplots.pdf")
# OutFile_path<-paste0(args[7],"/",OutFile_name)
# 
# write a PDF of the plot
pdf(file = paste0(args[7],"/alignment_barplots.pdf"),width = 8,height = 8)
layout(Raw_Reads_Matrix_matrix)
par(mar=c(0,0,4,0))
plot(1,type='n',axes=FALSE,xlab="",ylab="",main = "Sequencing Reads",cex.main=2)
legend("bottom",legend=c("Deduplicated","Duplicated","Unaligned"),fill=BARPLOT_COLORS,bty = "n",ncol=length(BARPLOT_COLORS),cex=1.0)
par(mar=c(6,8,0,3))
barplot(Dup_Raw_Reads_Matrix,horiz = T,col=BARPLOT_COLORS,border=NA,las=1,cex.names=0.7,xlab="Number of reads (millions)")
barplot(Dup_Pcnt_Reads_Matrix,horiz = T,col=BARPLOT_COLORS,border=NA,las=1,cex.names=0.7,xlab="Percent of reads")
dev.off()

# write a CSV of the final table
# # peel off the rownames into a separate vector
SampleName<-row.names(ReadsTable)
write.csv(x = cbind(SampleName,ReadsTable), file = paste0(args[7],"/alignment_stats_extended.csv"),quote = F,row.names = F)
