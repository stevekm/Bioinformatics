#!/bin/bash
set -x # verbose debugging output to stderr; turn this off when compiling report

## USAGE: peaks-per-sample.sh
# this script will tabulate the number of peaks per sample in each pipeline branch

# get the arguments
pipeline_peaks_dir="$1" # dir with the peak calling results
pipeline_annontate_peaks_dir="$2" # dir with the peak HOMER annotation results
main_outdir="$3" # outdir for the results

# load R version
module switch r/3.3.0

mkdir -p "$main_outdir"
cd "$main_outdir"

# make an array to hold the peak branches
declare -a index_array

# make an associative array to hold the paths to the files we will make
declare -A table_array

# get the indexes from the branch names for the peaks
index_array=($(find ${pipeline_peaks_dir}/results -maxdepth 1 -mindepth 1 -type d ! -name ".db" -exec basename {} \; | sort | tr '\n' ' '))

# make the aggregate barplots off the venn.txt
for key in ${!index_array[@]}; do
  echo "key is ${key}"
  
  # get the key value ; the branch
  branch=${index_array[${key}]}
  echo "branch is $branch"
  
  # make a subdir per branch
  tmp_outdir="${main_outdir}/${branch}"
  mkdir -p "$tmp_outdir"
  
  # set a file to hold the output peak stats table
  tmp_tablefile="${tmp_outdir}/peaks_stats.tsv"
  # print header
  echo -e "Peaks\tSample" > $tmp_tablefile
  
  # find all the peaks.bed files for the branch, 
  tmp_files=$(find "${pipeline_peaks_dir}" -path "*/${branch}/*" -name "peaks.bed" | tr '\n' ' ') # -exec wc -l {} \; # | tr ' ' '\n' # -print0  | xargs -0 $tmp_script
  
  # print the number of peaks and name for each sample into the table
  for i in $tmp_files; do
    tmp_sampleID="${i}"
    echo "tmp_sampleID is $tmp_sampleID"
    
    # get the sample ID from the dirname
    tmp_sampleID=$(basename $(dirname $tmp_sampleID))
    echo "new tmp_sampleID is $tmp_sampleID"
    
    # get the number of peaks 
    num_peaks=$(cat $i | wc -l)
    echo "num_peaks is $num_peaks"
    
    # print it to the table
    echo -e "${num_peaks}\t${tmp_sampleID}" >> "$tmp_tablefile"
    
  done
  
  # add the table file path to the tables array; use the branch name as the index
  table_array[$branch]="$(readlink -f $tmp_tablefile)"
  
  # pass the table file to R for plotting
  Rscript --slave --no-save --no-restore - "$tmp_tablefile" "$branch" <<EOF
  #!/usr/bin/env Rscript
  # R code for plotting the peaks table
  
  # ~~~~~ GET SCRIPT ARGS ~~~~~~~ #
  args <- commandArgs(TRUE); cat("Script args are:\n"); args
  
  # path to the peaks table file
  peaks_table_file <- args[1]
  
  # get the sample branch
  # peaks_branch <- basename(dirname(peaks_table_file))
  peaks_branch <- args[2]
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
  
  # load the file into a dataframe
  peaks_table_df<-read.table(peaks_table_file,header = TRUE,sep = "\t",stringsAsFactors = FALSE,check.names = FALSE)
  
  # convert the Sample column entries into rownames
  rownames(peaks_table_df) <- peaks_table_df[["Sample"]]
  
  
  # plot layout setup
  plot_layout_matrix<-structure(c(1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 
                                2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 
                                2L, 2L, 2L), .Dim = c(8L, 4L), .Dimnames = list(NULL, c("V1", 
                                                                                        "V2", "V3", "V4")))
  pdf(file = paste0(dirname(peaks_table_file),"/peaks_barplots.pdf"),width = 8,height = 8)
  # setup the panel layout
  layout(plot_layout_matrix)
  # need to set this for some reason
  # plot margins # c(bottom, left, top, right) # default is c(5, 4, 4, 2) + 0.1
  par(mar=c(0,0,5,0))
  # call blank plot to fill the first panel
  plot(1,type='n',axes=FALSE,xlab="",ylab="",main = peaks_branch,cex.main=1.5) 
  # set up the Legend in the first panel
  # legend("bottom",legend=colnames(overlap_df),bty = "n",cex=1.0) # fill=BARPLOT_COLORS,,ncol=length(BARPLOT_COLORS)
  # set some plot margin parameters to fit the names
  par(mar=c(5,14,0,2)+ 0.1) 
  barplot(t(peaks_table_df),
        # main=peaks_branch,
        cex.names = 0.7,
        horiz = T,
        # col=BARPLOT_COLORS,
        border=NA,
        las=1,
        # cex.names=Names_scale,
        xlab="Number of peaks",
        space=0.6
  ) 
  dev.off()
  save.image(file="data.Rdata",compress = TRUE)
EOF
  
done

# didn't need this part:
# # get all the tables we made
# for key in ${!table_array[@]}; do
#   echo "branch is ${key}"
#   
#   # get the key value ; the branch
#   tmp_table=${table_array[${key}]}
#   echo "tmp_table is $tmp_table"

# done

# THE END !!
exit
