#!/usr/bin/Rscript
# USAGE: motif_analysis_HOMER2_LocTable.R /path/to/known_motif_locations.txt
# READ IN THE HOMER Known Motif location table, parse it, and clean it up
# also output a mini version for testing; remove this later


# get commands passed to the script
args <- commandArgs()
# /path/to/stats.csv will be args[6]
# /path/to/HOMER_genes_list.txt will be args[7] 

##################

# return the args
Motif_Location_Table<-args[6]

# ~/projects//motif_analysis/Nkx2_2-1/known_motif_locations.txt

# read in the table
tmp<-read.table(file=Motif_Location_Table, sep = "\t",fill = T,quote="",header=T)
# tmp<-read.table(file="~/projects//motif_analysis/Nkx2_2-1/known_motif_locations.txt",sep = "\t",fill = T,quote="",header=T)

# fix the first column name
colnames(tmp)<-c("PeakID",colnames(tmp)[2:ncol(tmp)])

# fix the empty cells
tmp<-as.data.frame(apply(tmp, 2, function(x) gsub("^$|^ $", ".", x)))

# write out the objects
write.table(file=paste0(Motif_Location_Table,".tsv"),x = tmp,sep = "\t",na = "NA",row.names = F,quote = F)
write.table(file=paste0(Motif_Location_Table,"_mini.tsv"),x = tmp[1:100,],sep = "\t",na = "NA",row.names = F,quote = F)
