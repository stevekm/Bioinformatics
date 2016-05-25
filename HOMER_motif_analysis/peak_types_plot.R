#!/usr/bin/Rscript

## USAGE: Rscript peak_types_plot.R $outdir $peaks_type_table
## DESCRIPTION: create barplots to visualize the types of peaks, along with LaTeX formatted table
## 

# get the script arguments
args <- commandArgs(TRUE)

cat("R: dput args:") # for development & troubleshooting
cat("\n")
cat(dput(args))
cat("\n\n")
# 

OutDir<-args[1]
PeakTypesTable_file<-args[2]

# read in the table from the file
PeakTypesTable<-read.table(file = PeakTypesTable_file,header = FALSE,sep = "\t")
# dput(PeakTypesTable)
# structure(list(V1 = c(658L, 2L, 10L), V2 = structure(1:3, .Label = c("intron ", 
#                                                                      "non-coding ", "promoter-TSS "), class = "factor")), .Names = c("V1", 
#                                                                                                                                      "V2"), class = "data.frame", row.names = c(NA, -3L))
# > PeakTypesTable
# V1            V2
# 1 658       intron 
# 2   2   non-coding 
# 3  10 promoter-TSS 

# format the table
rownames(PeakTypesTable)<-PeakTypesTable$V2
PeakTypesTable<-PeakTypesTable["V1"]
colnames(PeakTypesTable)<-"NumberPeaks"
# > dput(PeakTypesTable)
# structure(list(NumberPeaks = c(658L, 2L, 10L)), .Names = "NumberPeaks", row.names = c("intron ", 
#                                                                                       "non-coding ", "promoter-TSS "), class = "data.frame")
# > PeakTypesTable
# NumberPeaks
# intron                658
# non-coding              2
# promoter-TSS           10

#
##
###
# get the total peaks
TotalPeaks<-sum(PeakTypesTable[["NumberPeaks"]])

#
##
###
# Calculate the percent 
PeakTypesTable[["PercentPeaks"]]<-signif((PeakTypesTable[["NumberPeaks"]]/TotalPeaks)*100,digits=3)

#
##
###
# save the values to be plotted into a transposed matrix, since thats what the barplot() likes
Peaks_Raw_Matrix<-t(as.matrix(PeakTypesTable["NumberPeaks"]))
# > dput(Peaks_Raw_Matrix)
# structure(c(658L, 2L, 10L), .Dim = c(1L, 3L), .Dimnames = list(
#   "NumberPeaks", c("intron ", "non-coding ", "promoter-TSS "
#   )))
# > Peaks_Raw_Matrix
# intron  non-coding  promoter-TSS 
# NumberPeaks     658           2            10
Peaks_Pcnt_Matrix<-t(as.matrix(PeakTypesTable["PercentPeaks"])) # don't actually need this..


#
##
###
# make a barplot
pdf(file = paste0(OutDir,"/peaks_type_barplots.pdf"),width = 8,height = 8)
barplot(Peaks_Raw_Matrix,main="Number of peaks",horiz = TRUE,las=1,cex.names=0.5)
dev.off()



#
##
###
# write a CSV of the final table
# # peel off the rownames into a separate vector
Types<-row.names(PeakTypesTable)
write.csv(x = cbind(Types,PeakTypesTable), file = paste0(OutDir,"/peaks_type_table.csv"),quote = FALSE,row.names = FALSE)

