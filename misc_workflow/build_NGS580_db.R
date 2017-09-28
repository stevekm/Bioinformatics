#!/usr/bin/env Rscript
# Script to create a SQLite database containing aggregated data for the NGS580 clinical database

library("data.table")
library("RSQLite")
library("digest")
library("reshape2")

# ~~~~~ PARAMS & SETTINGS ~~~~~ #
# a text file with the paths to the results of every clinical sequencing run analysis to use
results_dirs_paths_file <- "/ifs/data/molecpathlab/NGS580_WES-development/db/results_dirs_paths.txt"

# table of metadata from the NextSeq about runs; contains Matija's run ID's along with machine's ID's
NextSeq_run_index_file <- "/ifs/data/molecpathlab/quicksilver/run_index/NextSeq_run_index.csv"

# table with all NGS580 target genomic regions with ANNOVAR annotations
regions_annotations_file <- "/ifs/data/molecpathlab/NGS580_WES-development/db/target_regions.hg19_multianno.txt"

# .bed format file with fasta sequence for each region appended on each row
target_regions_fasta_file <- "/ifs/data/molecpathlab/NGS580_WES-development/db/target_regions_fasta.bed"

# list of genes included in the IonTorrent 50 gene panel
IonTorrent_reporter_panel_genes_file <- "/ifs/data/molecpathlab/NGS580_WES-development/db/IonTorrent_reporter_panel_genes.txt"

# database file to use
sqlite_file <- "/ifs/data/molecpathlab/NGS580_WES-development/db/NGS580.sqlite"



# ~~~~~ FUNCTIONS ~~~~~ #
msprintf <- function(fmt, ...) {
    message(sprintf(fmt, ...))
}

tsprintf <- function(fmt, ...){
    # print a formatted message with timestamp
    m <- sprintf(fmt, ...)
    message(sprintf('[%s] %s', format(Sys.time(), "%H:%M:%S"), m))
}

remove_empty_str <- function(x){
    # remove empty strings from character vector
    x <- x[which(! x %in% "")]
    return(x)
}

get_numlines <- function(input_file, skip = NA) {
    # count the number of lines in a file
    # skip = integer number to subtract from line count e.g. to skip header (doesn't actually prevent lines from being read)
    num_lines <- length(readLines(input_file))
    if(!is.na(skip)) num_lines <- num_lines - as.numeric(skip)
    return(num_lines)
}

chrom_rownames2cols <- function(df){
    # split rownames into separate columns for chromosome coordinates
    # chr10:100026989-100027328
    df_chrom <- as.data.frame(do.call(rbind, strsplit(rownames(df), ':')))
    df_chrom <- cbind(df_chrom[1], as.data.frame(do.call(rbind, strsplit(as.character(df_chrom$V2), '-'))))
    colnames(df_chrom) <- c("chrom", "start", "stop")
    df <- cbind(df_chrom, df)
    return(df)
}

sanitize_SQLite_colnames <- function(df){
    # clean the column names for use in SQLite
    bad_chars <- c('.', '-')
    for(bad_char in bad_chars){
        colnames(df) <- gsub(pattern = bad_char, replacement = '_', x = colnames(df), fixed = TRUE)
    }
    return(df)
}

add_uid <- function(df){
    # add a unique ID column 'uid' to a dataframe
    df[["uid"]] <- apply(X = df, MARGIN = 1, digest, algo = "xxhash64")
    return(df)
}

make_results_list <- function(paths){
    # make a list of runs & results IDs from paths to results dirs
    tsprintf("Making results list from run paths, searching for files associated with each run...")

        # empty list to hold results
    results_list <- list()
    
    # iterate over all the files passed
    for(i in seq_along(paths)){
        path <- paths[i]
        
        # parse the path
        run <- basename(dirname(path))
        results_ID <- basename(path)
        
        # find some files
        # standard pipeline output
        LoFreq_annot_all_file <- dir(path = path, pattern = "VCF-LoFreq-annot.all.txt", full.names = TRUE)
        GATK_HC_annot_all_file <- dir(path = path, pattern = "VCF-GATK-HC-annot.all.txt", full.names = TRUE)
        summary_combined_wes_file <- dir(path = path, pattern = "summary-combined.wes.csv", full.names = TRUE)
        
        # these files might not exist
        # use system 'find' due to recursive symlinks
        average_coverage_per_sample_file <- system(command = sprintf('find "%s" -name "*_average_coverage_per_sample.tsv"', path), intern = TRUE)
        
        
        
        
        
        # make sure file was found!
        if(length(average_coverage_per_sample_file) < 1){
            message(sprintf("WARNING: 'average_coverage_per_sample.tsv' not found for run %s, run will be skipped", path))
        } else {
            # add to the list
            tsprintf("Adding run %s to the results_list", path)
            results_list[[i]] <- list(run = run,
                                      results_ID = results_ID,
                                      path = path,
                                      LoFreq_annot_all_file = LoFreq_annot_all_file,
                                      GATK_HC_annot_all_file = GATK_HC_annot_all_file,
                                      average_coverage_per_sample_file = average_coverage_per_sample_file,
                                      summary_combined_wes_file = summary_combined_wes_file)
        }
    }
    return(results_list)
}

get_all_analysis_filepaths_table <- function(results_list){
    # search all the analysis directories and make a table with paths to every file per type
    tsprintf("Getting paths to all files for all analyses")
    
    data_list <- list()
    
    for(i in seq_along(results_list)){
        result <- results_list[[i]]
        run <- result[["run"]]
        results_ID <- result[["results_ID"]]
        
        path <- result[["path"]]
        
        analysis_dirs <- list.dirs(path, recursive = FALSE)
        
        files_list <- list()
        
        for(q in seq_along(analysis_dirs)){
            analysis_dir <- analysis_dirs[q]
            analysis_dir_name <- basename(analysis_dir)
            analysis_dir_files <- dir(analysis_dir, full.names = TRUE)
            df <- data.frame(path = analysis_dir_files)
            if(nrow(df) > 0){
                df[["dir"]] <- analysis_dir_name
                df[["file"]] <- basename(analysis_dir_files)
                df[["run"]] <- run
                df[["results_ID"]] <- results_ID
                
                files_list[[q]] <- df
            }
        }
        files_df <- as.data.frame(data.table::rbindlist(files_list))
        data_list[[i]] <- files_df
    }
    
    tsprintf("Building file paths dataframe...")
    return(as.data.frame(data.table::rbindlist(data_list)))
}

create_run_results_index <- function(results_list){
    # create an index of all run IDs, results IDs, and paths
    # return a df
    tsprintf("Making run_results_index")
    run_results_index <- data.frame()
    for(i in seq_along(results_list)){
        result <- results_list[[i]]
        
        
        run_df <- data.frame(run = result[["run"]],
                             results_ID = result[["results_ID"]],
                             path = result[["path"]])
        
        if(nrow(run_results_index) < 1){
            run_results_index <- run_df
        } else {
            run_results_index <- rbind(run_results_index, run_df)
        }
    }
    
    # remove dupes
    run_results_index <- run_results_index[! duplicated(run_results_index), ]
    
    return(run_results_index)
}

read_annotations <- function(file){
    # read the annotation files output by sns pipeline
    tsprintf("Reading file: %s", file)
    df <- read.delim(file = file, header = TRUE, sep = '\t', stringsAsFactors = FALSE, na.strings = ".")
    setnames(x = df, old = 'X.MUT', new = 'MUT')
    return(df)
}

read_all_annotations <- function(results_list, annot_file_list_index){
    # read all of the annotations files from the results_list into a single df
    # annot_file_list_index is the list index to read for each results entry
    tsprintf("Loading all annotations of type '%s' from all runs...", annot_file_list_index)
    # get vector of files to read in
    # annot_files <- sapply(X = results_list, "[[", annot_file_list_index)
    
    # empty list to hold results
    annot_dfs <- list()
    
    # iterate over the files
    for(i in seq_along(results_list)){
        result <- results_list[[i]]
        run <- result[["run"]]
        results_ID <- result[["results_ID"]]
        annot_file <- result[[annot_file_list_index]]
        
        # read the file into a dataframe
        df <- read_annotations(annot_file)

        # add extra columns
        df[["run"]] <- run
        df[["results_ID"]] <- results_ID

        annot_dfs[[i]] <- df
    }
    
    tsprintf("Building '%s' dataframe...", annot_file_list_index)
    return(as.data.frame(data.table::rbindlist(annot_dfs)))
}

add_all_annotations_to_db <- function(db_con, table_name, results_list, annot_file_list_index){
    tsprintf("Loading all annotations of type '%s' from all runs...", annot_file_list_index)
    # iterate over the files
    for(i in seq_along(results_list)){
        result <- results_list[[i]]
        run <- result[["run"]]
        results_ID <- result[["results_ID"]]
        annot_file <- result[[annot_file_list_index]]
        
        # read the file into a dataframe
        df <- read_annotations(annot_file)
        
        # add extra columns
        df[["run"]] <- run
        df[["results_ID"]] <- results_ID
        
        add_df_to_db(db_con = db_con,
                     table_name = table_name,
                     df = df)
    }
}

read_all_summary_combined <- function(results_list){
    # get all the summary-combined.wes.csv files and read them into a single table
    tsprintf("Loading all summary_combined files from all runs...")
    
    data_list <- list()
    
    for(i in seq_along(results_list)){
        result <- results_list[[i]]
        run <- result[["run"]]
        results_ID <- result[["results_ID"]]
        summary_combined_file <- result[["summary_combined_wes_file"]]

        df <- read.delim(file = summary_combined_file, header = TRUE, sep = ',', na.strings = 'X', check.names = TRUE, stringsAsFactors = FALSE)
        
        # fix colnames
        setnames(x = df, old = c("X.SAMPLE",
                                 "TRIM.SURVIVING.READS..",
                                 "MAPPED.READS..MQ10.",
                                 "MAPPED..",
                                 "CHIMERIC..",
                                 "DUPLICATES..",
                                 "ON.TARGET",
                                 "ON.TARGET.100BP.PAD",
                                 "ON.TARGET.500BP.PAD",
                                 "X._bases_above_10",
                                 "X._bases_above_50",
                                 "X._bases_above_100",
                                 "X._bases_above_500"),
                 new = c("SAMPLE",
                         "TRIM.SURVIVING.READS.PCNT",
                         "MAPPED.READS.MQ10",
                         "MAPPED.PCNT",
                         "CHIMERIC.PCNT",
                         "DUPLICATES.PCNT",
                         "ON.TARGET.PCNT",
                         "ON.TARGET.100BP.PAD.PCNT",
                         "ON.TARGET.500BP.PAD.PCNT",
                         "PCNT_bases_above_10",
                         "PCNT_bases_above_50",
                         "PCNT_bases_above_100",
                         "PCNT_bases_above_500"))
        colnames(df) <- gsub(pattern = '.', replacement = '_', x = colnames(df), fixed = TRUE)
        
        # add extra columns
        df[["run"]] <- run
        df[["results_ID"]] <- results_ID
        
        # generate unique IDs for each entry
        df <- add_uid(df)
        
        data_list[[i]] <- df
    }
    
    return(as.data.frame(data.table::rbindlist(data_list)))
}

read_all_avg_coverages <- function(results_list){
    # read in all the '_average_coverage_per_sample.tsv' files to create a single dataframe with all data from all runs
    
    tsprintf("Reading in all average coverage files from all runs")
    # cov_files <- sapply(X = results_list, "[[", "average_coverage_per_sample_file")
    
    # empty list to hold values
    cov_dfs <- list()
    
    # column names for non-samples
    id_cols <- c("chrom", "start", "stop", "region")
    
    # for(cov_file in cov_files){
    for(i in seq_along(results_list)){
        result <- results_list[[i]]
        run <- result[["run"]]
        results_ID <- result[["results_ID"]]
        cov_file <- result[["average_coverage_per_sample_file"]]
        
        tsprintf("Reading file: %s", cov_file)
        
        df <- read.delim(file = cov_file, header = TRUE, sep = '\t', row.names = 1, check.names = FALSE)
        
        # add region columns and chrom cols
        df[["region"]] <- rownames(df)
        df <- chrom_rownames2cols(df)
        
        # melt the df to long format
        df <- reshape2::melt(df,
                             id.vars = id_cols,
                             variable.name = "sample",
                             value.name = "coverage")
        
        # add run info
        df[["run"]] <- run
        df[["results_ID"]] <- results_ID
        
        cov_dfs[[i]] <- df
    }
    
    tsprintf("Building coverages dataframe...")
    return(as.data.frame(data.table::rbindlist(cov_dfs)))
}


add_all_coverages_to_db <- function(db_con, table_name = 'Coverage_per_region', results_list){
    tsprintf("Adding all average coverage files from all runs to database")
    
    
    # column names for non-samples
    id_cols <- c("chrom", "start", "stop", "region")
    
    # for(cov_file in cov_files){
    for(i in seq_along(results_list)){
        result <- results_list[[i]]
        run <- result[["run"]]
        results_ID <- result[["results_ID"]]
        cov_file <- result[["average_coverage_per_sample_file"]]
        
        tsprintf("Reading file: %s", cov_file)
        
        df <- read.delim(file = cov_file, header = TRUE, sep = '\t', row.names = 1, check.names = FALSE)
        
        # add region columns and chrom cols
        df[["region"]] <- rownames(df)
        df <- chrom_rownames2cols(df)
        
        # melt the df to long format
        df <- reshape2::melt(df,
                             id.vars = id_cols,
                             variable.name = "sample",
                             value.name = "coverage")
        
        # add run info
        df[["run"]] <- run
        df[["results_ID"]] <- results_ID
        
        add_df_to_db(db_con = db_con,
                     table_name = table_name,
                     df = df)
    }
    
}

read_IT50_genes <- function(path){
    # read in the list of genes for the IonTorrent 50 gene panel
    # need to do extra processing on the file
    genes <- readLines(path)
    
    # remove blank entries
    genes <- remove_empty_str(genes)
    
    # add missing gene
    genes <- c(genes, "GNA11")
    
    # remove bad gene
    genes <- genes[which(! genes == "GNA 11")]
    return(genes)
}

create_regions_annotations_fasta_table <- function(target_regions_fasta, regions_annotations, IonTorrent_genes){
    # combine the annotated regions with the fasta sequences,
    # do some calculations,
    # add info about IT50 gene panel entries
    
    # clean up the df's
    colnames(target_regions_fasta) <- c("region", "fasta")
    regions <- do.call(paste, c(regions_annotations[c("Chr", "Start")], sep = ":"))
    regions <- paste(regions, regions_annotations[["End"]], sep = "-")
    regions_annotations[["region"]] <- regions
    regions_annotations_fasta <- merge(x = regions_annotations, y = target_regions_fasta, by = "region")
    
    # calculate region stats
    regions_annotations_fasta[["fasta"]] <- as.character(regions_annotations_fasta[["fasta"]])
    regions_annotations_fasta[["fasta"]] <- toupper(regions_annotations_fasta[["fasta"]])
    
    regions_annotations_fasta[["G_total"]] <- sapply(regmatches(regions_annotations_fasta[["fasta"]],
                                                                gregexpr("G",
                                                                         regions_annotations_fasta[["fasta"]])), length)
    regions_annotations_fasta[["C_total"]] <- sapply(regmatches(regions_annotations_fasta[["fasta"]],
                                                                gregexpr("C",
                                                                         regions_annotations_fasta[["fasta"]])), length)
    regions_annotations_fasta[["T_total"]] <- sapply(regmatches(regions_annotations_fasta[["fasta"]],
                                                                gregexpr("T",
                                                                         regions_annotations_fasta[["fasta"]])), length)
    regions_annotations_fasta[["A_total"]] <- sapply(regmatches(regions_annotations_fasta[["fasta"]],
                                                                gregexpr("A",
                                                                         regions_annotations_fasta[["fasta"]])), length)
    
    # count the total length o the sequence
    regions_annotations_fasta[["fasta_total"]] <- sapply(regions_annotations_fasta[["fasta"]], nchar)
    
    regions_annotations_fasta[["GC_total"]] <- regions_annotations_fasta[["G_total"]] + regions_annotations_fasta[["C_total"]]
    regions_annotations_fasta[["GC_content"]] <- regions_annotations_fasta[["GC_total"]] / regions_annotations_fasta[["fasta_total"]]
    
    
    regions_annotations_fasta[["in_IT50"]] <- FALSE
    regions_annotations_fasta[which(regions_annotations_fasta[["Gene.refGene"]] %in% IonTorrent_genes), "in_IT50"] <- TRUE
    
    return(regions_annotations_fasta)
}


make_regions_df <- function(target_regions_fasta_file, regions_annotations_file, IonTorrent_reporter_panel_genes_file){
    df <- create_regions_annotations_fasta_table(
        # load target fastas
        target_regions_fasta = read.delim(file = target_regions_fasta_file,
                                          header = FALSE,
                                          sep = '\t',
                                          na.strings = '.'),
        # load annotated target regions
        regions_annotations = read.delim(file = regions_annotations_file,
                                         header = TRUE,
                                         sep = '\t',
                                         stringsAsFactors = FALSE,
                                         na.strings = '.'),
        # load the IonTorrent 50 gene panel genes
        IonTorrent_genes = read_IT50_genes(path = IonTorrent_reporter_panel_genes_file)
    )
    return(df)
}


load_NextSeq_run_index <- function(NextSeq_run_index_file){
    NextSeq_run_index <- read.delim(file = NextSeq_run_index_file, header = TRUE, sep = ',')
    return(NextSeq_run_index)
}

# SQLite Functions
df_sqlite_names <- function(df){
    # creates character string of the column names for use in SQLite statement
    sql_names <- sprintf('"%s"', names(df))
    return(sql_names)
}

df_sqlite_param_names <- function(df){
    # creates character string of the column names for use in SQLite prepared statement
    sql_names <- sprintf('@%s', names(df))
    return(sql_names)
}


create_table_with_primary_key <- function(db_con, table_name, df, key_colname = "uid"){
    # create a table in a SQLite database and set the primary key for the table
    tsprintf("Creating table '%s' in database and adding all entries...", table_name)
    if(! table_name %in% dbListTables(db_con)){
        sql_names <- df_sqlite_names(df)
        sql <- sprintf("CREATE TABLE %s(%s, primary key(%s))",
                       table_name,
                       paste(sql_names, collapse = ", "),
                       key_colname)
        # message(sql)
        dbGetQuery(db_con, sql)
        dbWriteTable(db_con, table_name, df, append = TRUE, row.names = FALSE)
    } else {
        tsprintf("Table '%s' already exists in database, skipping...", table_name)
    }
}


insert_or_ignore <- function(db_con, table_name, df, speedup = FALSE){
    # insert entries from the df into the db if they're not already present
    tsprintf("Inserting new entries into table '%s' in database...", table_name)
    sql_names <- df_sqlite_names(df)
    param_names <- df_sqlite_param_names(df)
    
    sql <- sprintf('INSERT OR IGNORE INTO %s(%s) VALUES (%s)',
                   table_name,
                   paste(sql_names, collapse = ", "),
                   paste(param_names, collapse = ", "))
    # message(sql)
    if(speedup == TRUE){
        dbGetQuery(db_con, "PRAGMA synchronous = OFF") 
        dbGetQuery(db_con, "PRAGMA journal_mode = OFF")
    }
    
    dbBegin(db_con)
    dbGetPreparedQuery(db_con, sql, bind.data=df)
    dbCommit(db_con)
    # 1: RSQLite::dbGetPreparedQuery() is deprecated, please switch to DBI::dbGetQuery(params = bind.data).
}

sqlite_count_table_rows <- function(db_con, table_name){
    sql <- sprintf('SELECT COUNT(*) FROM %s', table_name)
    row_count <- unlist(dbGetQuery(db_con, sql))[["COUNT(*)"]]
    return(row_count)
}


add_df_to_db <- function(db_con, table_name, df){
    # add dataframe to SQLite database with primary key set
    # if the table does not exist; create it
    # if the table exists; insert only new entries
    
    tsprintf("Adding '%s' to database...", table_name)
    
    # add 'uid' unique ID column and clean the column names
    tsprintf("Preparing dataframe for insertion into database...")
    df <- add_uid(df)
    df <- sanitize_SQLite_colnames(df)

    if( ! table_name %in% dbListTables(db_con)){
        # add table if its not present
        tsprintf("Table '%s' does not exist, adding table and data...", table_name)
        create_table_with_primary_key(db_con = db_con,
                                      table_name = table_name,
                                      df = df)
    } else {
        # if the table already exists, add its entries
        num_rows_start <- sqlite_count_table_rows(db_con = db_con, table_name = table_name)
        
        tsprintf("Table '%s' already exists in database with %s rows", table_name, num_rows_start)
        
        tsprintf("Starting insert operation for dataframe containing %s rows...", nrow(df))

        insert_or_ignore(db_con = db_con, table_name = table_name, df = df)
        
        num_rows_end <- sqlite_count_table_rows(db_con = db_con, table_name = table_name)
        
        tsprintf("Finished insertion for table '%s', table now contains %s rows", table_name, num_rows_end)
    }
}






# ~~~~ RUN ~~~~~ #
# load the NextSeq run index
NextSeq_run_index <- load_NextSeq_run_index(NextSeq_run_index_file)

# load the paths to results dirs
results_dirs_paths <- readLines(results_dirs_paths_file, skipNul = TRUE)
results_dirs_paths <- remove_empty_str(results_dirs_paths)

# make list of runs
results_list <- make_results_list(paths = results_dirs_paths)


# open db connection
mydb <- dbConnect(RSQLite::SQLite(), sqlite_file)

# add items to database
# small datasets
add_df_to_db(db_con = mydb,
             table_name = 'analysis_files',
             df = get_all_analysis_filepaths_table(results_list = results_list))

add_df_to_db(db_con = mydb,
             table_name = 'run_results_index',
             df = create_run_results_index(results_list))

add_df_to_db(db_con = mydb,
             table_name = 'summary_combined_wes',
             df = read_all_summary_combined(results_list))

add_df_to_db(db_con = mydb,
             table_name = 'regions_annotations_fasta',
             df = make_regions_df(target_regions_fasta_file,
                                  regions_annotations_file,
                                  IonTorrent_reporter_panel_genes_file))

# problems here with duplicate entries or something
# add_df_to_db(db_con = mydb,
#              table_name = 'NextSeq_run_index',
#              df = NextSeq_run_index)


# very large datasets
add_all_annotations_to_db(db_con = mydb, results_list = results_list, annot_file_list_index = "LoFreq_annot_all_file", table_name = 'LoFreq_annotations')

add_all_annotations_to_db(db_con = mydb, results_list = results_list, annot_file_list_index = "GATK_HC_annot_all_file", table_name = 'GATK_HC_annotations')

add_all_coverages_to_db(db_con = mydb, results_list = results_list)


dbDisconnect(mydb)




