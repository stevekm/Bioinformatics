#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
vcf_tsv_file <- args[1]
annovar_txt_file <- args[2]
avinput_file <- args[3]
output_file <- args[4]
# example:
# vcf_tsv_file <- "NC-HAPMAP.updated.tsv"
# annovar_txt_file <- "NC-HAPMAP.hg19_multianno.txt"
# avinput_file <- 'NC-HAPMAP.avinput.reformat.tsv'
# output_file <- "output.tsv"


# ~~~~~~ CUSTOM FUNCTIONS ~~~~~~~ #
read.ANNOVAR.vcf_txt <- function(file, na_string = '.', drop_Otherinfo = TRUE){
    # read the '_multianno.txt' file produced by ANNOVAR table_annovar.pl using the '--vcfinput' arg
    # tab-delimited text file has more columns than column names
    # output dataframe with colnames present
    
    # read in the file without headers; headers in first row
    tmp_annovar_df <- read.delim(file = file, 
                                 header = FALSE, 
                                 sep = '\t', 
                                 check.names = FALSE, 
                                 fill = TRUE, stringsAsFactors = FALSE, 
                                 na.strings = na_string)
    # get the colnames that are present
    av_colnames <- tmp_annovar_df[1,][which(! is.na(tmp_annovar_df[1,]) & ! tmp_annovar_df[1,] == "" )]
    # create new df without colnames
    tmp_annovar_df2 <- tmp_annovar_df[2:nrow(tmp_annovar_df),]
    # add the colnames
    names(tmp_annovar_df2)[1:length(av_colnames)] <- av_colnames
    
    # drop columns after 'Otherinfo', if its present
    if(isTRUE(drop_Otherinfo)){
        if (any(grepl(pattern = 'Otherinfo', x = names(tmp_annovar_df2)))){
            Otherinfo_index <- which(names(tmp_annovar_df2) == 'Otherinfo')
            tmp_annovar_df2 <- tmp_annovar_df2[, names(tmp_annovar_df2)[1:Otherinfo_index - 1]]
        }
    }
    
    return(tmp_annovar_df2)
}

all_equal <- function(...){
    # returns TRUE or FALSE if all the elements passed are equal
    arguments <- list(...)
    TF <- vapply(1:(length(arguments)-1),
                 function(n) identical(arguments[[n]], arguments[[n+1]]),
                 logical(1))
    return(all(TF))
}

# ~~~~~~ RUN ~~~~~~ #
annovar <- read.ANNOVAR.vcf_txt(file = annovar_txt_file)
vcf_table <- read.delim(file = vcf_tsv_file, header = TRUE, sep = '\t')
avinput <- read.delim(file = avinput_file, header = TRUE, sep = '\t')
save.image("loaded.Rdata")

# check that all df's have the same number of rows
if(! all_equal(nrow(annovar), nrow(avinput), nrow(vcf_table)) ) {
    print("ERROR: Loaded dataframes have unequal numbers of rows")
    quit(status = 1)
}

# merge tables
merged_df <- Reduce(function(x, y){ merge(x, y, all = TRUE) }, list(avinput, vcf_table, annovar))
save.image("merged.Rdata")

# check that all df's have the same number of rows
if(! all_equal(nrow(annovar), nrow(avinput), nrow(vcf_table), nrow(merged_df)) ) {
    print("ERROR: Dataframes have unequal numbers of rows after merging")
    quit(status = 1)
}

write.table(x = merged_df, file = output_file, sep = '\t', quote = FALSE, na = '.', row.names = FALSE, col.names = TRUE)
save.image("finish.Rdata")