#!/usr/bin/Rscript
# USAGE: annovar2redcap.R /path/to/annotation/file.csv /path/to/original_file.vcf ANNOVAR_version /path/to/CONVERTED_TSV OUTPUT_GENE_TABLE_BASENAME.csv RUN_NAME /path/to/REDCap_data_dictionary.csv
# This script will convert the ANNOVAR annotated CSV files into a format that REDCap wants
# args[6]: /path/to/annotation/file.csv
# args[7]: /path/to/original_file.vcf # get the basename..
# args[8]: ANNOVAR_version
# args[9]: VCF converted to TSV by VCFLIB <https://github.com/ekg/vcflib>

# get commands passed to the script
args <- commandArgs()

# self test.. print some stuff
print("R: input is: ")
print(args[6:length(args)])

# # for testing purposes:
# print("ARGS dput results here:")
# dput(args)
# # # # # 
# 
# args<-c( copy some args here from the dput output)
# quit()



# get the ANNOVAR table previously created
ANNOVAR_table<-read.csv(file=args[6])
FileName_original<-basename(args[7])
ANNOVAR_version<-args[8]
Converted_TSV<-args[9]
Converted_TSV_df<-read.table(file = Converted_TSV,header = T,sep = "\t")
OUTPUT_GENE_TABLE_BASENAME<-args[10]
RUN_NAME<-args[11]
REDCap_Data_Dictionary<-read.csv(file = args[12])

# setup the table function
ANNOVAR_geneTable<-function(ANNOVAR_table, # the ANNOVAR output
                            Converted_TSV_df, # needed for the allele frequency
                            REDCap_Data_Dictionary){ # REDCap database setup file that has already been created for the project
  
  #
  # # 
  # # #
  # if no rows, exit
  #
  if (nrow(ANNOVAR_table) == 0) { print("R: EMPTY ANNOVAR TABLE, EXITING")
    quit() }
  if (nrow(Converted_TSV_df) == 0) { print("R: EMPTY VCF TSV TABLE, EXITING")
    quit() }
  
  # 
  # # 
  # # # 
  # Set up the ANNOVAR_table for latter
  #
  # # Get allele frequency from the Converted_TSV_df and add that ANNOVAR_table
  # # concatenate the chr, Pos/Start, Ref, Alt columns into one entry
  # # for both df's, then merge both to get the AF from TSV into the ANNOVAR_table
  # # make sure that the AF column exists in the TSV table; allele frequency
  if("AF" %in% colnames(Converted_TSV_df)) # if there is a column labeled AF
  {
    # create tmp ID to merge the df's on 
    Converted_TSV_df$tmp_join_id<-with(Converted_TSV_df,paste(CHROM,POS,REF,ALT,sep = "-"))
    
    # keep only the tmp ID and the af's
    Converted_TSV_df<-Converted_TSV_df[,c("tmp_join_id","AF")]
    
    # create tmp ID to merge the df's on 
    ANNOVAR_table$tmp_join_id<-with(ANNOVAR_table,paste(Chr,Start,Ref,Alt,sep = "-"))
    
    # merge the two tables
    ANNOVAR_table<-merge(ANNOVAR_table,Converted_TSV_df,by = "tmp_join_id",all.x = TRUE)
    
    # rename this column
    colnames(ANNOVAR_table)[colnames(ANNOVAR_table) == "AF" ]<-"allele_frequency"
    
    
  } else { # fill the AF col with NA's if it doesn't exist
    
    allele_frequency<-rep("NA",times=nrow(ANNOVAR_table))
    ANNOVAR_table<-cbind(ANNOVAR_table,allele_frequency)
  }
  
  
  # # set up the ANNOVAR_table columns for the gene table later
  # 
  # # adjust ANNOVAR_table colnames to match REDCap criteria
  # # getthe correct columns as per the Data_Dictionary
  
  
  # # # REDCap criteria for colnames:
  # Variables/field names must consist of ONLY lower-case letters, numbers, and underscores.
  # # so replace . with _ , lowercase all, colnames
  colnames(ANNOVAR_table)<-gsub(pattern = '.',replacement = "_",x = colnames(ANNOVAR_table),fixed = T)
  colnames(ANNOVAR_table)<-tolower(x = colnames(ANNOVAR_table))
  
  # get the list of columns for the Gene_table from the REDCap Data Dictionary
  Desired_Columns<-as.character(REDCap_Data_Dictionary$Variable...Field.Name)
  
  # check to make sure that each column we want is in the ANNOVAR table; if not, create it and fill with "."
  for(i in Desired_Columns[!Desired_Columns %in% colnames(ANNOVAR_table)]){ # all the desrd. cols. not in the table
    # make new col in table, fill with '.'
    ANNOVAR_table[[paste(i)]]<-rep(".",times=nrow(ANNOVAR_table))
  }
  
  
  
  #
  # # 
  # # # 
  # Make the Gene Table
  # 
  
  # select the columns we want, make new df
  Gene_table<-ANNOVAR_table[,paste(Desired_Columns)]
  
  # VARIANT FILTERING
  # # # DISABLE THE FILTERING.. for now
  # remove unknown and synonymous SNV entries
  # Gene_table<-Gene_table[!(Gene_table$ExonicFunc.refGene %in% c("unknown", "synonymous SNV")),]
  # keep only the exon based entries
  # Gene_table<-Gene_table[(Gene_table$Func.refGene %in% c("exonic","exonic;splicing")),]
  
  # Sort by the scores # TURN THIS OFF TOO
  # Gene_table<-Gene_table[ with(Gene_table,order(-(as.numeric(SIFT_score)),-(as.numeric(Polyphen2_HDIV_score)))), ]
  # Sort by 
  
  # REDCap wants a CSV; but there are entries in the data frame table that have commas... replace with " "
  Gene_table<-as.data.frame(apply(Gene_table,MARGIN=c(1,2),function(x){
    gsub(",", " ", x,fixed = T)
  }))
  
  # fill in the empty columns
  Gene_table$sample_name<-rep(FileName_original,times=nrow(Gene_table))
  Gene_table$run_name<-rep(RUN_NAME,times=nrow(Gene_table))
  Gene_table$unique_id<-paste(RUN_NAME,FileName_original,seq(1:nrow(Gene_table)),sep=":")
  Gene_table$annovar_version<-rep(ANNOVAR_version,times=nrow(Gene_table))
  Gene_table$gene_table_file<-rep(OUTPUT_GENE_TABLE_BASENAME,times=nrow(Gene_table))
  Gene_table$annovar_complete<-rep(2,times=nrow(Gene_table))
  
  
  # clean up the clinvar entry a little..
  Gene_table$clinvar_20150629<-gsub(x = Gene_table$clinvar_20150629,pattern = ";",replacement = "; ")
  return(Gene_table)
}

genetable_hashR<-function(Gene_table){
  #
  # # 
  # in development; lets make hash values for each variant
  # # do both variant ID's independent and dependent on sample/run
  library(digest)
  tmp_Gene_table2<-Gene_table
  # generate a hash based just on the variant information; potentially useful to match same variants between samples
  tmp_Gene_table2$variant_hash<-apply(tmp_Gene_table2[,c("chr","start","end","ref","alt")],1,digest)
  
  # generate a hash based on the variant info + sample info
  tmp_Gene_table2$unique_hash<-apply(tmp_Gene_table2[,c("run_name","sample_name","chr","start","end","ref","alt")],1,digest)
  # not sure what to use this for..
  return(tmp_Gene_table2)
}






# run the function on the input files
Gene_table<-ANNOVAR_geneTable(ANNOVAR_table = ANNOVAR_table,
                              Converted_TSV_df = Converted_TSV_df,
                              REDCap_Data_Dictionary = REDCap_Data_Dictionary)
# write the output to a file
Out_File<-paste0(dirname(args[6]),"/",gsub(pattern = ".hg19_multianno.csv",x = basename(args[6]),replacement = "_"),"gene_table.csv")
print("R: writing file:")
print(Out_File)
write.csv(Gene_table,file=Out_File,quote = F,row.names = F)


# experimental gene table with hash values..
Gene_table<-genetable_hashR(Gene_table = Gene_table)
# write the output to a file
Out_File<-paste0(dirname(args[6]),"/",gsub(pattern = ".hg19_multianno.csv",x = basename(args[6]),replacement = "_"),"gene_table_hashR.csv")
print("R: writing file:")
print(Out_File)
write.csv(Gene_table,file=Out_File,quote = F,row.names = F)
