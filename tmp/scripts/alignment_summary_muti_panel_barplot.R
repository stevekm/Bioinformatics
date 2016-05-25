# This code creates a multi-panel plot with barplots displaying the total and percentage values for alignment reads
##################
# Set up the data
# sample data frame with reads
ReadsTable<-structure(list(Total.Reads = c(43743601.2, 37888316.7, 44664123.66, 39804329.76, 36708833.46, 49666748.82, 47457348.9, 48866464.92), 
                           Aligned.Reads = c(33114730.74, 28309936.98, 34478776.02, 30168289.32, 27943990.62, 33808528.26, 33249888.36, 37155596.4), 
                           De.duplicated.alignments = c(32004797.76, 26982520.5, 29728794.42, 28707722.94, 26717704.26, 25559734.98, 22385229.36, 34029148.92)), 
                      .Names = c("Total.Reads", "Aligned.Reads", "De.duplicated.alignments"), 
                      row.names = c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6", "Sample7", "Sample8"), 
                      class = "data.frame")
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


#####################
# Set up the plots

BARPLOT_COLORS<-c("blue","purple","red")

# setup the matrix for the plot layout
Raw_Reads_Matrix_matrix<-structure(c(1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 1L, 2L, 
                                     2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 
                                     3L, 3L, 3L, 3L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L),
                                   .Dim = c(11L,4L), 
                                   .Dimnames = list(NULL, c("V1", "V2", "V3", "V4")))
# looks like this
# > Raw_Reads_Matrix_matrix
# V1 V2 V3 V4
# [1,]  1  1  1  1
# [2,]  2  2  2  2
# [3,]  2  2  2  2
# [4,]  2  2  2  2
# [5,]  2  2  2  2
# [6,]  2  2  2  2
# [7,]  3  3  3  3
# [8,]  3  3  3  3
# [9,]  3  3  3  3
# [10,]  3  3  3  3
# [11,]  3  3  3  3

# setup the panel layout
layout(Raw_Reads_Matrix_matrix)
# layout.show(max(Raw_Reads_Matrix_matrix))  # use this to preview the layout

par(mar=c(0,0,4,0)) # need to set this for some reason
# on mar: A numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot. 
# The default is c(5, 4, 4, 2) + 0.1.

# call blank plot to fill the first panel
plot(1,type='n',axes=FALSE,xlab="",ylab="",main = "Sequencing Reads",cex.main=2) 

# set up the Legend in the first panel
legend("bottom",legend=c("Deduplicated","Duplicated","Unaligned"),fill=BARPLOT_COLORS,bty = "n",ncol=length(BARPLOT_COLORS),cex=1.0)

# set the mar for the remaining panels
par(mar=c(6,5,0,4))
# create barplot for the two matrices
barplot(Dup_Raw_Reads_Matrix,horiz = T,col=BARPLOT_COLORS,border=NA,las=1,cex.names=1,xlab="Number of reads (millions)")
barplot(Dup_Pcnt_Reads_Matrix,horiz = T,col=BARPLOT_COLORS,border=NA,las=1,cex.names=1,xlab="Percent of reads")

# write a PDF version
pdf("alignment_summary_stats_reads_barplot.pdf",width = 8,height = 8)
layout(Raw_Reads_Matrix_matrix)
par(mar=c(0,0,4,0))
plot(1,type='n',axes=FALSE,xlab="",ylab="",main = "Sequencing Reads",cex.main=2)
legend("bottom",legend=c("Deduplicated","Duplicated","Unaligned"),fill=BARPLOT_COLORS,bty = "n",ncol=length(BARPLOT_COLORS),cex=1.0)
par(mar=c(6,5,0,4))
barplot(Dup_Raw_Reads_Matrix,horiz = T,col=BARPLOT_COLORS,border=NA,las=1,cex.names=1,xlab="Number of reads (millions)")
barplot(Dup_Pcnt_Reads_Matrix,horiz = T,col=BARPLOT_COLORS,border=NA,las=1,cex.names=1,xlab="Percent of reads")
dev.off()
