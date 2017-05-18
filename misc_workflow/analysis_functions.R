# functions for the workflow
read_diffbind_sheet <- function(input_file, verbose = FALSE) {
    # read in a DiffBind sheet, return the df
    if(verbose == TRUE) message(sprintf("Reading in file:\n%s", input_file))
    df <- read.table(file = input_file, header = TRUE, sep = ",")
    return(df)
}

# 
# write_BED_file <- function(df, output_file, verbose = FALSE){
#     # output the first 3 columns of the dataframe as a BED formatted file; chrom start stop
#     if(verbose == TRUE) message(sprintf("Writing output to file:\n%s", output_file))
#     write.table(x = df[1:3], file = output_file, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
# }

write_BED_file <- function(df, output_file, verbose = TRUE){
    # output the first 3 columns of the dataframe as a BED formatted file; chrom start stop
    if(verbose == TRUE) message(sprintf("Writing output to file:\n%s", output_file))
    # make the parent dir if doesn't exist
    if(! dir.exists(dirname(output_file))) dir.create(path = dirname(output_file), recursive = TRUE)
    write.table(x = df[1:3], file = output_file, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
}


remove_duplicated_rows <- function(df, verbose = FALSE){
    # remove duplicated rows in the dataframe
    start_length <- nrow(df)
    df <- df[! duplicated(df), ]
    end_length <- nrow(df)
    if(verbose == TRUE) message(sprintf("Number of duplicate rows removed:\n%s", start_length - end_length))
    return(df)
}

make_output_filename <- function(name, suffix){
    return(sprintf("%s_%s", name, suffix))
}

get_sample_ID_from_filename <- function(filename){
    # "/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2017-12-31/project_notes/integrated_analysis/source_data/diffbind_sheets/H3K9AC.D-vs-R.blocking.p100.csv"
    # output: H3K9AC
    sample_ID <- gsub(pattern = "^(.*?)\\..*$", replacement = "\\1", perl = TRUE, x = basename(filename))
    return(sample_ID)
}

save_DiffBind_loci <- function(diffbind_file, output_dir){
    # diffbind_file is a character vector element
    # output_dir is character vector path to outdir
    df <- read_diffbind_sheet(diffbind_file, verbose = TRUE)
    df <- remove_duplicated_rows(df[1:3], verbose = TRUE)
    sample_ID <- get_sample_ID_from_filename(diffbind_file)
    write_BED_file(df = df, verbose = TRUE,
                   output_file = file.path(output_dir,
                                           make_output_filename(sample_ID, "diffbind_regions.bed")))
}

read_microarray_sheet <- function(input_file){
    df <- read.delim(file = input_file, header = TRUE, sep = '\t')
    rownames(df) <- df[[1]]
    df <- df[2:ncol(df)]
    return(df)
}

rename_microarray_df_cols <- function(df){
    colnames(df) <- gsub(pattern = '.Exp', replacement = '', x = colnames(df))
    return(df)
}

write_microarray_table <- function(df, output_file){
    write.table(x = df, quote = FALSE, sep = '\t', row.names = TRUE, col.names=NA, file = output_file)
}

read_methylation_sheet <- function(input_file){
    df <- read.delim(file = input_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
    rownames(df) <- df[[1]]
    # df <- df[2:ncol(df)]
    return(df)
}

split_df_col <- function(df, colname_to_split, split_char, new_colnames){
    # split a column in a dataframe based on a character, return the df with the new split col
    split_col_df <- as.data.frame(do.call(rbind, strsplit(as.character(df[[colname_to_split]]), split_char)))
    colnames(split_col_df) <- new_colnames
    return(cbind(split_col_df, df[, colnames(df)[! colnames(df) %in% colname_to_split] ]))
}


clean_colname <- function(old_name){
    return(gsub(pattern = '-', replacement = '.', x = as.character(old_name)))
}

replace_colname <- function(df, old_colname, new_colname, clean_colnames = FALSE){
    if(clean_colnames == TRUE){
        old_colname <- clean_colname(old_colname)
        new_colname <- clean_colname(new_colname)
    }
    colnames(df)[grep(pattern = old_colname, x = colnames(df))] <- new_colname
    return(df)
}

read_diffbind_methlyation_regions_sheet <- function(input_file){
    df <- read.delim(file = input_file, header = FALSE, sep = '\t')
    head(df)
    colnames(df) <- c("diffbind.chrom", "diffbind.start", "diffbind.stop", "methylation.chrom", "methylation.start", "methylation.stop")
    return(df)
    
}

read_BED_file_with_geneIDs <- function(input_file){
    df <- read.delim(file = input_file, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
    colnames(df) <- c("chrom", "start", "stop", "gene")
    return(df)
}


read_methylation_regions_sheet <- function(input_file){
    df <- read_BED_file_with_geneIDs(input_file)
    df <- remove_duplicated_rows(df)
    return(df)
}

make_diffbind_data_regions_list <- function(diffbind_region_files, diffbind_data_files){
    histone_marks <- c("H3K27AC", "H3K27ME3", "H3K4ME3", "H3K9AC", "H3K9ME3")
    diffbind_data_list <- list()
    for(mark in histone_marks){
        region_file <- grep(pattern = mark, x = diffbind_region_files, value = TRUE)
        data_file <- grep(pattern = mark, x = diffbind_data_files, value = TRUE)
        
        diffbind_data_list[[mark]] <- list( regions = read_diffbind_methlyation_regions_sheet(region_file),
                                            data = read_diffbind_sheet(data_file))
    }
    return(diffbind_data_list)
}

read_methylation_sheet_fix_chromLocs <- function(input_file){
    df <- read_methylation_sheet(input_file)
    df <- split_df_col(df, colname_to_split = "X", split_char = ":", new_colnames = c("chrom", "loc"))
    df <- split_df_col(df, colname_to_split = "loc", split_char = "-", new_colnames = c("start", "stop"))
    df[["chrom"]] <- as.character(df[["chrom"]])
    return(df)
}

no_chrom_df <- function(df){
    return(df[ colnames(df)[! colnames(df) %in% c("start", "stop", "chrom")] ])
}

make_character_cols <- function(df, cols_to_convert){
    for(col in cols_to_convert){
        df[[col]] <- as.character(df[[col]])
    }
    return(df)
}

quantile_normalize_df <- function(df){
    # quantile normalization of a dataframe
    df <- setnames(x = as.data.frame(normalize.quantiles(as.matrix(df))), colnames(df))
    return(df)
}
