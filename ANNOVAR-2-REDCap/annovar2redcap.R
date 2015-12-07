#!/usr/bin/Rscript
# USAGE: annovar2redcap.R /path/to/annotation/file.csv /path/to/original_file.vcf ANNOVAR_version /path/to/CONVERTED_TSV OUTPUT_GENE_TABLE_BASENAME.csv RUN_NAME
# args[6]: /path/to/annotation/file.csv
# args[7]: /path/to/original_file.vcf # get the basename..
# args[8]: ANNOVAR_version
# args[9]: VCF converted to TSV by VCFLIB <https://github.com/ekg/vcflib>

# get commands passed to the script
args <- commandArgs()

# self test.. print some stuff
print("R: input is: ")
print(args[6:length(args)])

# get the ANNOVAR table previously created
ANNOVAR_table<-read.csv(file=args[6])
FileName_original<-basename(args[7])
ANNOVAR_version<-args[8]
Converted_TSV<-args[9]
Converted_TSV_df<-read.table(file = Converted_TSV,header = T,sep = "\t")
OUTPUT_GENE_TABLE_BASENAME<-args[10]
RUN_NAME<-args[11]

# if no rows, exit
if (nrow(ANNOVAR_table) == 0) { print("R: EMPTY ANNOVAR TABLE, EXITING")
  quit() }
if (nrow(Converted_TSV_df) == 0) { print("R: EMPTY VCF TSV TABLE, EXITING")
  quit() }

# make sure that the AF column exists in the TSV table (allele frequency)
if("AF" %in% colnames(Converted_TSV_df))
{
  # # create tmp ID to merge the df's on
  Converted_TSV_df$tmp_join_id<-with(Converted_TSV_df,paste(CHROM,POS,REF,ALT,sep = "-"))
  Converted_TSV_df<-Converted_TSV_df[,c("tmp_join_id","AF")]
  
  ANNOVAR_table$tmp_join_id<-with(ANNOVAR_table,paste(Chr,Start,Ref,Alt,sep = "-"))
  
  # # merge the two tables
  ANNOVAR_table<-merge(ANNOVAR_table,Converted_TSV_df,by = "tmp_join_id",all.x = TRUE)
  
  # rename this column
  colnames(ANNOVAR_table)[colnames(ANNOVAR_table) == "AF" ]<-"allele_frequency"
  
} else {
  # fill the AF col with NA's
  allele_frequency<-rep("NA",times=nrow(ANNOVAR_table))
  ANNOVAR_table<-cbind(ANNOVAR_table,allele_frequency)
}

# set a list of columns that we want; include allele freq
Desired_Columns<-c("Gene.refGene","Chr","Start","End","Ref","Alt","Func.refGene","ExonicFunc.refGene",
                   "AAChange.refGene","allele_frequency","SIFT_score","SIFT_pred","Polyphen2_HDIV_score","Polyphen2_HDIV_pred","cosmic68")

# check to make sure that each column we want is in the ANNOVAR table; if not, create it and fill with "."
for(i in Desired_Columns[!Desired_Columns %in% colnames(ANNOVAR_table)]){ # all the desrd. cols. not in the table
  # make new col in table, fill with '.'
  ANNOVAR_table[[paste(i)]]<-rep(".",times=nrow(ANNOVAR_table))
}

# select the columns we want, make new df
Gene_table<-ANNOVAR_table[paste(Desired_Columns)]

# # # DISABLE THE FILTERING.. for now
# remove unknown and synonymous SNV entries
# Gene_table<-Gene_table[!(Gene_table$ExonicFunc.refGene %in% c("unknown", "synonymous SNV")),]
# keep only the exon based entries
# Gene_table<-Gene_table[(Gene_table$Func.refGene %in% c("exonic","exonic;splicing")),]

# Sort by the scores # TURN THIS OFF TOO
# Gene_table<-Gene_table[ with(Gene_table,order(-(as.numeric(SIFT_score)),-(as.numeric(Polyphen2_HDIV_score)))), ]
# Sort by 

# for REDCap, ; Variables/field names must consist of ONLY lower-case letters, numbers, and underscores.
# so replace . with _ , lowercase all, colnames
colnames(Gene_table)<-gsub(pattern = '.',replacement = "_",x = colnames(Gene_table),fixed = T)
colnames(Gene_table)<-tolower(x = colnames(Gene_table))

# REDCap wants a CSV; but there are entries in the data frame table that have commas... replace with " "
Gene_table<-as.data.frame(apply(Gene_table,MARGIN=c(1,2),function(x){
  gsub(",", " ", x,fixed = T)
}))

# need a unique ID for REDCap first field...
#  add filename/samplename column
gene_table_id<-seq(1:nrow(Gene_table))
sample_name<-rep(FileName_original,times=nrow(Gene_table))
run_name<-rep(RUN_NAME,times=nrow(Gene_table))
# concatenate the sample names and the ID numbers
unique_id<-paste(run_name,sample_name,gene_table_id,sep=":")

Gene_table<-cbind(sample_name,Gene_table)
Gene_table<-cbind(run_name,Gene_table)
Gene_table<-cbind(unique_id,Gene_table)


# add column for ANNOVAR version
annovar_version<-rep(ANNOVAR_version,times=nrow(Gene_table))
Gene_table<-cbind(Gene_table,annovar_version)

# add a column for the final gene_table.csv file
gene_table_file<-rep(OUTPUT_GENE_TABLE_BASENAME,times=nrow(Gene_table))
Gene_table<-cbind(Gene_table,gene_table_file)

# add column for ANNOVAR entry completion, since REDCap wants this; use value '2'
annovar_complete<-rep(2,times=nrow(Gene_table))
Gene_table<-cbind(Gene_table,annovar_complete)

Out_File<-paste0(dirname(args[6]),"/",gsub(pattern = ".hg19_multianno.csv",x = basename(args[6]),replacement = "_"),"gene_table.csv")
print("R: writing file:")
print(Out_File)
write.csv(Gene_table,file=Out_File,quote = F,row.names = F)

