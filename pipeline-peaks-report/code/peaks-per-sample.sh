#!/bin/bash
set -x # verbose debugging output to stderr; turn this off when compiling report

## USAGE: peaks-per-sample.sh
# this script will tabulate the number of peaks per sample in each pipeline branch

# load R version
module switch r/3.3.0
# dir with the peak calling results
pipeline_peaks_dir="/ifs/data/smithlab/SmithLab-ChIPSeq_2016-12-31/pipeline/peaks/results"
pipeline_annontate_peaks_dir="/ifs/data/smithlab/SmithLab-ChIPSeq_2016-12-31/pipeline/homer_annotatepeaks/results"
pipeline_alignstats_sheet="/ifs/data/smithlab/SmithLab-ChIPSeq_2016-12-31/pipeline/align-stats/results/align-stats.standard/align.by_sample.bowtie2/all-samples/alignment_stats_extended.csv"

# outdir for all overlap sets
# /ifs/data/smithlab/SmithLab-ChIPSeq_2016-06-06/project_notes/peaks-per-sample_report/analysis_pipeline
# main_outdir="/ifs/data/smithlab/SmithLab-ChIPSeq_2016-06-06/project_notes/peaks-per-sample"
main_outdir="/ifs/home/kellys04/projects/SmithLab-ChIPSeq_2016-12-31/analysis_dir/project_notes/peaks-per-sample_report_3/output"
mkdir -p "$main_outdir"
cd "$main_outdir"

# make an array to hold the peak branches
declare -a index_array

# make an associative array to hold the paths to the files we will make
declare -A table_array

# get the indexes from the branch names for the peaks
index_array=($(find ${pipeline_peaks_dir}/ -maxdepth 1 -mindepth 1 -type d ! -name ".db" -exec basename {} \; | sort | tr '\n' ' '))

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
  Rscript --slave --no-save --no-restore - "$tmp_tablefile" "$branch" "$pipeline_alignstats_sheet" <<EOF
  #!/usr/bin/env Rscript
  # R code for plotting the peaks table
  
  # ~~~~~ GET SCRIPT ARGS ~~~~~~~ #
  args <- commandArgs(TRUE); cat("Script args are:\n"); args
  
  # path to the peaks table file
  peaks_table_file <- args[1]
  
  # get the sample branch
  # peaks_branch <- basename(dirname(peaks_table_file))
  peaks_branch <- args[2]
  
  align_stats_file <- args[3]
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
  
  
  # ~~~~~ LOAD PEAKS FILE ~~~~~~~ #
  # load the peaks file into a dataframe
  peaks_table_df<-read.table(peaks_table_file,header = TRUE,sep = "\t",stringsAsFactors = FALSE,check.names = FALSE)
  
  # convert the Sample column entries into rownames
  rownames(peaks_table_df) <- peaks_table_df[["Sample"]]
  
  # re order based on rownames
  peaks_table_df <- peaks_table_df[order(rownames(peaks_table_df)),]
  
  
  # ~~~~~ LOAD ALIGN FILE ~~~~~~~ #
  align_stats_df <- read.delim(align_stats_file, header = TRUE, sep=",",row.names = 1)
  align_stats_df
  
  # re order the rownames
  align_stats_df <- align_stats_df[order(rownames(align_stats_df)),]
  
  # save the values to be plotted into a transposed matrix, since thats what the barplot() likes
  # # first just get the columns we want
  Dup_Raw_Reads_df<-align_stats_df[,which(colnames(align_stats_df) %in% c("De.duplicated.alignments","Duplicated","Unaligned.Reads")) ] 
  # reorder the columns because R is dumb
  Dup_Raw_Reads_df<-Dup_Raw_Reads_df[c("De.duplicated.alignments","Duplicated","Unaligned.Reads")]
  
  Dup_Raw_Reads_Matrix<-t(as.matrix(Dup_Raw_Reads_df))
  # # divid the number of reads by 1million
  Dup_Raw_Reads_Matrix<-signif(Dup_Raw_Reads_Matrix/1000000,digits = 4)
  
  # # first just get the columns we want
  Dup_Pcnt_Reads_df<-align_stats_df[,which(colnames(align_stats_df) %in% c("Percent.De.dup.Reads","Percent.Dup","Pcnt.Unaligned.Reads")) ]
  # reorder
  Dup_Pcnt_Reads_df<-Dup_Pcnt_Reads_df[c("Percent.De.dup.Reads","Percent.Dup","Pcnt.Unaligned.Reads")]
  Dup_Pcnt_Reads_Matrix<-t(as.matrix(Dup_Pcnt_Reads_df))
  
  
  # make align stats plot
  # Set up the plots
  BARPLOT_COLORS<-c("blue","purple","red")
  # setup the matrix for the plot layout
  Raw_Reads_Matrix_matrix<-structure(c(1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 1L, 2L, 
                                       2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 
                                       3L, 3L, 3L, 3L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L),
                                     .Dim = c(11L,4L), 
                                     .Dimnames = list(NULL, c("V1", "V2", "V3", "V4")))
  Raw_Reads_Matrix_matrix <- Raw_Reads_Matrix_matrix[(rowSums(Raw_Reads_Matrix_matrix < 3) > 0), , drop = FALSE]
  
  # calculate some values for the figure margins
  # based on the nchar() of the longest rowname, divided by a value
  mar_divisor<-2.0 # # smaller number = larger margin ; might have to adjust this manually per-project
  mar_widthLeft<-signif(
    max(4.1,
        max(nchar(row.names(align_stats_df)))/mar_divisor),
    4) # 65 is too much # this may not work in RStudio but should work in pdf() or command line R 
  
  
  # calculate value for the plot label scaling factor
  # Names_scale<-min(0.7, 
  #                max(nchar(row.names(align_stats_df)))*.0075) # works alright up to 88 char's samplenames # this doesn't work as well, needs more tweaking
  # ^ this also causes tiny labels for short names, too small, need both max and min cutoffs??
  Names_scale<-0.7
  # scaling factor for space between bars
  Space_scale<-max(0.2, # default setting
                   nrow(align_stats_df)*.01) # needs to work with up to 61
  
  
  cat("The longest sample name is ",max(nchar(row.names(align_stats_df))),sep = "\n")
  cat("The number of samples is ",nrow(align_stats_df),sep = "\n")
  
  
  cat("mar_widthLeft is ",mar_widthLeft,"",sep = "\n")
  cat("Names_scale is ",Names_scale,"",sep = "\n")
  cat("Space_scale is ",Space_scale,sep = "\n")
  
  # write a PDF of the plot
  # pdf(file = paste0(OutDir,"/alignment_barplots",mar_divisor,"-",mar_widthLeft,".pdf"),width = 8,height = 8) # ORIGINAL
  pdf(file = paste0(dirname(peaks_table_file),"/alignment_barplots.pdf"),width = 8,height = 8)
  
  # setup the panel layout
  layout(Raw_Reads_Matrix_matrix) 
  # need to set this for some reason
  par(mar=c(0,0,4,0))
  # call blank plot to fill the first panel
  plot(1,type='n',axes=FALSE,xlab="",ylab="",main = "Sequencing Reads",cex.main=2) 
  # set up the Legend in the first panel
  legend("bottom",legend=c("Deduplicated","Duplicated","Unaligned"),fill=BARPLOT_COLORS,bty = "n",ncol=length(BARPLOT_COLORS),cex=1.0)
  # plot margins # c(bottom, left, top, right) # default is c(5, 4, 4, 2) + 0.1
  # par(mar=c(6,max(4.1,max(nchar(row.names(align_stats_df)))/1.5),0,3)+ 0.1) # ORIGINAL 
  par(mar=c(6,mar_widthLeft,0,3)+ 0.1) 
  
  # create barplot for the two matrices
  # barplot(Dup_Raw_Reads_Matrix,horiz = T,col=BARPLOT_COLORS,border=NA,las=1,cex.names=0.7,xlab="Number of reads (millions)") 
  # barplot(Dup_Pcnt_Reads_Matrix,horiz = T,col=BARPLOT_COLORS,border=NA,las=1,cex.names=0.7,xlab="Percent of reads")
  barplot(Dup_Raw_Reads_Matrix,horiz = T,col=BARPLOT_COLORS,border=NA,las=1,cex.names=Names_scale,xlab="Number of reads (millions)",space=Space_scale) 
  
  dev.off()
  
  
  
  
  
  
  # ~~~~~ START PEAKS PLOT ~~~~~~~ #
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
  par(mar=c(5,16,0,2)+ 0.1) 
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
  # p <- recordPlot()
  dev.off()
  
  
  
  
  # ~~~~~ START DUAL PLOT ~~~~~~~ #
  
  dual_plot_matrix <- structure(c(1, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 
  2, 1, 2, 2, 2, 2, 2, 3, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 3, 4, 
  4, 4, 4, 4, 3, 4, 4, 4, 4, 4), .Dim = c(6L, 8L), .Dimnames = list(
    NULL, c("V1", "V2", "V3", "V4", "", "", "", "")))
  dual_plot_matrix
  
  # start the dual plot
  pdf(file = paste0(dirname(peaks_table_file),"/dual_barplots.pdf"),width = 16,height = 8)
  # setup the panel layout
  layout(dual_plot_matrix) 
  
  # call blank plot to fill the first panel
  plot(1,type='n',axes=FALSE,xlab="",ylab="",main = "Sequencing Reads",cex.main=2) 
  # set up the Legend in the first panel
  par(mar=c(0,0,5,0))
  legend("bottom",legend=c("Deduplicated","Duplicated","Unaligned"),fill=BARPLOT_COLORS,bty = "n",ncol=length(BARPLOT_COLORS),cex=1.0)
  
  # plot margins # c(bottom, left, top, right) # default is c(5, 4, 4, 2) + 0.1
  # par(mar=c(6,mar_widthLeft,0,3)+ 0.1) 
  # par(mar=c(6,16,0,3)+ 0.1) 
  par(mar=c(5,16,0,2)+ 0.1)
  # create alignemnt barplot 
  barplot(Dup_Raw_Reads_Matrix,horiz = T,
        col=BARPLOT_COLORS,
        border=NA,
        las=1,
        cex.names=Names_scale,
        xlab="Number of reads (millions)"
        #,space=Space_scale
        ) 
  
  
  # start peaks plot
  par(mar=c(0,0,5,0))
  # call blank plot to fill the first panel
  plot(1,type='n',axes=FALSE,xlab="",ylab="",main = peaks_branch,cex.main=1.5) 
  # set up the Legend in the first panel
  
  par(mar=c(5,16,0,2)+ 0.1) 
  barplot(t(peaks_table_df),
        # main=peaks_branch,
        cex.names = 0.7,
        horiz = T,
        # col=BARPLOT_COLORS,
        border=NA,
        las=1,
        # cex.names=Names_scale,
        xlab="Number of peaks"
        #,space=0.6
  )
  dev.off()
  
  
  
  
  # ~~~~~ MERGE DATA SETS ~~~~~~~ #
  # make 'sampleID' columns in both data sets
  align_stats_df[["SampleID"]] <- rownames(align_stats_df)
  peaks_table_df[["SampleID"]] <- rownames(peaks_table_df)
  
  peaks_align_merge_df <- base::merge(align_stats_df,peaks_table_df,by="SampleID",all=TRUE)
  
  # only keep these columns
  peaks_align_merge_df <- peaks_align_merge_df[,c("SampleID","Total.reads","Aligned.reads","De.duplicated.alignments","Duplicated","Unaligned.Reads","Peaks")]
  
  # melt it into long format
  library("reshape2")
  peaks_align_merge_long_df <- reshape2::melt(peaks_align_merge_df,id.vars="SampleID",variable.name="variable",value.name="value")
  
  
  
  # ~~~~~ MAKE MERGED PLOTS ~~~~~~~ #
  library("grid")
  library("ggplot2")
  library("plyr")
  library("gridExtra")
  
  theme_set(theme_bw())
  # set up the sample ID's along the center
  g.mid<-ggplot(peaks_align_merge_long_df,aes(x=1,y=SampleID))+geom_text(aes(label=SampleID),size = 2)+
    # geom_segment(aes(x=0.94,xend=0.95,yend=SampleID))+ # geom_segment(aes(x=0.94,xend=0.96,yend=SampleID))+
    # geom_segment(aes(x=1.04,xend=1.05,yend=SampleID))+
    ggtitle("Samples")+
    ylab(NULL)+
    scale_x_continuous(expand=c(0,0),limits=c(0.94,1.065))+
    theme(axis.title=element_blank(),
          # panel.grid=element_line(size = 10),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.background=element_blank(),
          axis.text.x=element_text(color=NA),
          axis.ticks.x=element_line(color=NA),
          plot.margin = unit(c(1,-1,1,-1), "mm"))
  
  # subset for just the peaks values, order by rownames e.g. Sample IDs
  peaks_align_merge_peaksonly <- droplevels(subset(peaks_align_merge_long_df, variable=="Peaks"))
  # make sure they are ordered by sample ID
  # peaks_table_df <- peaks_table_df[order(rownames(peaks_table_df)),]
  peaks_align_merge_peaksonly <- peaks_align_merge_peaksonly[order(peaks_align_merge_peaksonly[["SampleID"]]),]
  
  # with(peaks_align_merge_peaksonly,order(-SampleID,value,))
  
  # make barplot for just the peaks
  g1 <- ggplot(data = peaks_align_merge_peaksonly, aes(x = SampleID, y = value)) +
    geom_bar(stat = "identity") + ggtitle("Peaks") +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          plot.margin = unit(c(1,-1,1,0), "mm")) +
    scale_y_reverse() + coord_flip()
  
  # subset for just the align stats reads values to plot:
  # c("Deduplicated","Duplicated","Unaligned")
  peaks_align_merge_readsonly <- peaks_align_merge_long_df[peaks_align_merge_long_df[["variable"]] %in% c("De.duplicated.alignments","Duplicated","Unaligned.Reads"),]
  
  # make barplot for just the align stats
  g2 <- ggplot(data = peaks_align_merge_readsonly, aes(x = SampleID, y = value, fill=variable)) +xlab(NULL)+
    geom_bar(stat = "identity") +  ggtitle("Total.reads") +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = unit(c(1,0,1,-1), "mm")
    ) + coord_flip() + scale_fill_manual(values=c("blue", "purple", "red"))
  # , legend.box = "horizontal", legend.position = "top") + labs(fill="Reads: ") + coord_flip()
  
  # put it all together
  gg1 <- ggplot_gtable(ggplot_build(g1))
  gg2 <- ggplot_gtable(ggplot_build(g2))
  gg.mid <- ggplot_gtable(ggplot_build(g.mid))
  
  
  # out_file_path <- "/ifs/home/kellys04/projects/SmithLab-ChIPSeq_2016-12-31/analysis_dir/project_notes/peaks-per-sample_report_3/output/peaks.by_group.macs_broad/gg_peaks_alignemtn_barplots.pdf"
  pdf(file = paste0(dirname(peaks_table_file),"/dual_ggbarplots.pdf"),width = 10,height = 8.5)
  grid.arrange(gg1,gg.mid,gg2,ncol=3,widths=c(0.3,0.2,0.5)) # widths=c(0.4,0.2,0.4)
  dev.off()

  
  # ~~~~~ SAVE PEAK ALIGN STATS ~~~~~~~ #
  # save the table to a TSV table, to print in the report
  peaks_align_merge_df_save <- peaks_align_merge_df[,c("SampleID","Peaks","De.duplicated.alignments")]
  write.table(x = peaks_align_merge_df_save,file = paste0(dirname(peaks_table_file),"/peak_align_stats.tsv"),quote = FALSE,sep = '\t',row.names = FALSE,col.names = TRUE)
  
  # save R session
  save.image(file=paste0(dirname(peaks_table_file),"/peak-stats.Rdata"),compress = TRUE)
  
  
  
  
  
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
exit
exit

# notes:
tmp_script="/ifs/data/smithlab/SmithLab-ChIPSeq_2016-06-06/project_notes/peaks-per-sample_report/code/peaks-per-sample.sh"
chmod +x $tmp_script
$tmp_script

