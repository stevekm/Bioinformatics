#!/usr/bin/Rscript
# USAGE: annovar2redcap.R path/to/annotation/file.csv
# this will be args[6]

# get commands passed to the script
args <- commandArgs()

# self test.. print some stuff
print(args[6:length(args)])
print(length(args))
print(basename(args[6]))
print(dirname(args[6]))
# 
# get the ANNOVAR table previously created
ANNOVAR_table<-read.csv(file=args[6])

# select the columns we want
Gene_table<-ANNOVAR_table[c("Gene.refGene","Chr","Start","End","Ref","Alt","Func.refGene","ExonicFunc.refGene",
                            "AAChange.refGene","SIFT_score","SIFT_pred","Polyphen2_HDIV_score","Polyphen2_HDIV_pred","cosmic68")]


# remove unknown and synonymous SNV entries
Gene_table<-Gene_table[!(Gene_table$ExonicFunc.refGene %in% c("unknown", "synonymous SNV")),]

# keep only the exon based entries
Gene_table<-Gene_table[(Gene_table$Func.refGene %in% c("exonic","exonic;splicing")),]

# Sort by the scores
Gene_table<-Gene_table[ with(Gene_table,order(-(as.numeric(SIFT_score)),-(as.numeric(Polyphen2_HDIV_score)))), ]

# for REDCap, ; Variables/field names must consist of ONLY lower-case letters, numbers, and underscores.
# so replace . with _ , lowercase all, colnames
colnames(Gene_table)<-gsub(pattern = '.',replacement = "_",x = colnames(Gene_table),fixed = T)
colnames(Gene_table)<-tolower(x = colnames(Gene_table))

# REDCap wants a CSV; but there are entries in the data frame table that have commas... replace with " "
Gene_table<-as.data.frame(apply(Gene_table,MARGIN=c(1,2),function(x){
  gsub(",", " ", x,fixed = T)
}))

# need a unique ID for REDCap first field... 
# gene_table_id<-paste0(Gene_table$chr,"_",Gene_table$start,"_",Gene_table$ref,"_",Gene_table$alt)
# gene_table_id<-gsub(pattern = " ",replacement = "",x = gene_table_id)
gene_table_id<-seq(1:nrow(Gene_table))
Gene_table<-cbind(gene_table_id,Gene_table)

Out_File<-paste0(dirname(args[6]),"/",gsub(pattern = ".hg19_multianno.csv",x = basename(args[6]),replacement = "_"),"gene_table.csv")
write.csv(Gene_table,file=Out_File,quote = F,row.names = F)
