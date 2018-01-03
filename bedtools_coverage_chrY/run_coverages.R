#!/usr/bin/env Rscript

# get the coverage of select regions from chrY in the exome sequencing targets
# run bedtools coverage and capture the output for all samples in an aggregate df
# http://bedtools.readthedocs.io/en/latest/content/tools/coverage.html
# example command:
# bedtools coverage -sorted -g "hg19.chrom.sizes" -a "IDT_chrY.bed" -b "BAM-BWA/15-1067.bam" | cut -f7-

# ~~~~~ SETUP ~~~~~ # 

# hg19 chrom sizes file
hg19_chrom_file <- "hg19.chrom.sizes"

# $ fetchChromSizes hg19 > hg19.chrom.sizes
# $ grep -v '_' hg19.chrom.sizes
# chr1	249250621
# chr2	243199373
# chr3	198022430
# chr4	191154276
# chr5	180915260
# chr6	171115067
# chr7	159138663
# chrX	155270560
# chr8	146364022
# chr9	141213431
# chr10	135534747
# chr11	135006516
# chr12	133851895
# chr13	115169878
# chr14	107349540
# chr15	102531392
# chr16	90354753
# chr17	81195210
# chr18	78077248
# chr20	63025520
# chrY	59373566
# chr19	59128983
# chr22	51304566
# chr21	48129895
# chrM	16571

# target intervals file to find coverages at
intervals_file <- "IDT_chrY.bed"
intervals_colnames <- c("chr", "start", "stop", "interval_ID", 'NA', "strand")

# get the paths to the .bam files to use for finding coverages
bam_files <- dir("BAM-BWA", pattern = '.bam$', full.names = TRUE)
# get their sample IDs from the filename
sampleIDs <- gsub(pattern = '.dd.ra.rc.bam', replacement = '', x = basename(bam_files))

sample_files <- setNames(object = bam_files, nm = sampleIDs)
#                               15-1067                                   248                                   249 
# "BAM-GATK-RA-RC/15-1067.dd.ra.rc.bam"     "BAM-GATK-RA-RC/248.dd.ra.rc.bam"     "BAM-GATK-RA-RC/249.dd.ra.rc.bam" 
#                                   250                                   252                                   253 
#     "BAM-GATK-RA-RC/250.dd.ra.rc.bam"     "BAM-GATK-RA-RC/252.dd.ra.rc.bam"     "BAM-GATK-RA-RC/253.dd.ra.rc.bam" 
#                                   254                                   255                                   256 
#     "BAM-GATK-RA-RC/254.dd.ra.rc.bam"     "BAM-GATK-RA-RC/255.dd.ra.rc.bam"     "BAM-GATK-RA-RC/256.dd.ra.rc.bam" 
#                                   257                                   258                                   259 
#     "BAM-GATK-RA-RC/257.dd.ra.rc.bam"     "BAM-GATK-RA-RC/258.dd.ra.rc.bam"     "BAM-GATK-RA-RC/259.dd.ra.rc.bam" 
#                                   260                                   261 
#     "BAM-GATK-RA-RC/260.dd.ra.rc.bam"     "BAM-GATK-RA-RC/261.dd.ra.rc.bam" 


# colnames of the extra columns of output added from bedtools coverage
bedtools_cov_colnames <- c("num_reads", "num_bases_covered", "interval_size", "fraction_bases_covered")





# ~~~~~ FUNCTIONS ~~~~~ # 
bedtools_cov_command <- function(a_file, b_file, module_cmd = 'module load bedtools/2.26.0', hg19_chrom_file = "hg19.chrom.sizes"){
    # build the system command for running bedtools coverage
    # a_file: file of regions to find depth for
    # b_file: file of coverages/reads/etc to count depth in
    # http://bedtools.readthedocs.io/en/latest/content/tools/coverage.html
    # example command:
    # bedtools coverage -sorted -g "hg19.chrom.sizes" -a "IDT_chrY.bed" -b "BAM-BWA/15-1067.bam" | cut -f7-
    
    bedtools_cov_cmd <- sprintf('%s; bedtools coverage -sorted -g "%s" -a "%s" -b "%s"', module_cmd, hg19_chrom_file, a_file, b_file)
    return(bedtools_cov_cmd)
}

char2df <- function(x, split_str = '\t'){
    # split a character vector on a string value, then rbind it into a dataframe
    # split on tabs by default
    df <- as.data.frame(do.call("rbind", strsplit(x, split_str, fixed = TRUE)))
    return(df)
}




# ~~~~~ RUN ~~~~~ #
# empty df to hold the coverages for all samples
intervals_cov_df <- data.frame()
for(i in seq_along(sample_files)){
    bam_file <- sample_files[i]
    sampleID <- names(sample_files)[i]
    
    # build the system command
    bedtools_cov_cmd <- bedtools_cov_command(a_file = intervals_file, b_file = bam_file)
    
    # run the system command, capture output in character vector
    cov_out <- system(bedtools_cov_cmd, intern = TRUE)
    
    # convert char vec to df by splitting on tabs
    cov_df <- char2df(x = cov_out)
    
    # add the colnames
    sample_bedtools_cols <- sprintf('%s_%s', sampleID, bedtools_cov_colnames)
    colnames(cov_df) <- c(intervals_colnames, sample_bedtools_cols)
    
    
    # check if the final df is empty
    if(ncol(intervals_cov_df) < 1){
        # add the base interval cols
        intervals_cov_df <- cov_df[, intervals_colnames]
    }
    
    # add the sample cols
    for(sample_bedtools_col in sample_bedtools_cols){
        intervals_cov_df[[sample_bedtools_col]] <- cov_df[[sample_bedtools_col]]
    }
}


# fix 'gene' col in the df; interval_ID col has gene ID between '()'
geneIDs <- gsub(pattern = '^.*\\(', replacement = '', x = intervals_cov_df[["interval_ID"]], perl = TRUE)
geneIDs <- gsub(pattern = '\\).*$', replacement = '', x = geneIDs, perl = TRUE)
intervals_cov_df[["gene"]] <- geneIDs

# save the output
write.table(x = intervals_cov_df, file = 'samples_intervals_num_reads_chrY_BAM-BWA.tsv', quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
saveRDS(object = intervals_cov_df, file = 'intervals_cov_df.RDS', ascii = TRUE)
