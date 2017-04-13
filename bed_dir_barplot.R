#!/usr/bin/env Rscript

## USAGE: bed_dir_barplot.R -n "Sample ID" -d /path/to/bedfiles/dir
## DESCRIPTION: This script will create a barplot for the number of entries (lines) in each BED file
## in the provided directory. 


# ~~~~~ LOAD PACKAGES ~~~~~ #
library("optparse")


# ~~~~~ CUSTOM FUNCTIONS ~~~~~ #
wc_l <- function(filepath){
    # count the number of lines in each file, using 'wc -l' because sometimes BED files might be too large for R
    return(as.numeric(system(sprintf("cat %s | wc -l", filepath), intern = TRUE)))
}


# ~~~~~ GET SCRIPT ARGS ~~~~~ #
option_list <- list(
    make_option(c("-n", "--name"), type = "character", dest = "name",
                help="Sample ID or Sample Name"),
    make_option(c("-x", "--exclude"), type = "character", default = NA, dest = "exclude",
                help="Filename pattern to exclude from plots"),
    make_option(c("-d", "--dir"), type = "character", dest = "dir",
                help="Directory containing BED files to plot")
)
opt <- parse_args(OptionParser(option_list=option_list))

sampleID <- opt$name
bed_dir <- opt$dir
exclude_pattern <- opt$exclude
# save.image()


# ~~~~~ GET BED FILES ~~~~~ #
bed_files <- sort(dir(bed_dir, pattern = "*.bed", full.names = TRUE)) 
if( ! is.na(exclude_pattern)){
    bed_files <- bed_files[grep(pattern = as.character(exclude_pattern), x = basename(bed_files), invert = TRUE)]
}

# ~~~~~ MAKE PLOT ~~~~~ #
pdf(file = file.path(bed_dir, "peak_total_numbers.pdf"), width = 8, height = 8)
par(mar=c(5,12,4,2)) # increase y-axis margin; default is  5.1 4.1 4.1 2.1
barplot(sapply(X = bed_files, FUN = wc_l), 
        names.arg = basename(bed_files), 
        main = sprintf("%s Total Peaks", sampleID), 
        horiz = TRUE, 
        cex.names = 0.7, 
        las=1)
dev.off()
