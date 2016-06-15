#!/usr/bin/env Rscript

## USAGE: peaks_overlap_aggregatoR.R /path/to/venn.txt /path/to/venn.txt /path/to/venn.txt ... ... 
## This script will process venn output from HOMER mergePeaks venn.txt for many samples and aggregate
## the data into a single table and barplot, saved into the current working dir
## only use this on pairwise venn's e.g. mergePeaks output on only 2 bed files

# example bash usage:
# tmp_script="/path/to/code/peaks_overlap_aggregatoR.R"
# peak_dir="/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2016-03-31/project_notes/peak_overlap/manuscript_overlaps_report/peak_overlaps"
# find "$peak_dir" -path "*/peaks.by_sample.macs_broad/*" -name "venn.txt" -print0 | xargs -0 $tmp_script


# ~~~~~ LOAD PACKAGES ~~~~~~~ #
# library('VennDiagram')
# library('gridExtra')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 


# ~~~~~ GET SCRIPT ARGS ~~~~~~~ #
args <- commandArgs(TRUE); cat("Script args are:\n"); args
# dput(args)
# get the plot filename from the first arg filepath
plot_filename<-basename(dirname(args[1]))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 


# ~~~~~ PARSE VENN TABLE FUNCTION ~~~~~~~ #
get_venn_entry <- function(venn_path) {
  # NOTE: the renaming commands here are customized for my specific use case; adjust accordingly as needed
  venn_table_file<-venn_path
  # read in the venn text
  venn_table_df<-read.table(venn_table_file,header = TRUE,sep = "\t",stringsAsFactors = FALSE,check.names = FALSE)
  # venn_table_df
  # looks like this:
#   BVI-D-H3K27AC.bed	BVI-R-H3K27AC.bed	Total	Name
#   X	41572	BVI-R-H3K27AC.bed
#   X		1011	BVI-D-H3K27AC.bed
#   X	X	13420	BVI-D-H3K27AC.bed|BVI-R-H3K27AC.bed
  
  # get the venn categories; remove Total and Name
  venn_categories<-colnames(venn_table_df)[!colnames(venn_table_df) %in% c("Total","Name")] 
  
  # venn_categories
  num_categories<-length(venn_categories)
  
  # make a summary table; strip the columns with X, etc.
  venn_summary<-venn_table_df[!colnames(venn_table_df) %in% venn_categories]
  # looks like this:
#   Total Name
#   28999 AGK-R-H3K27AC.bed
#   19518 AGK-D-H3K27AC.bed
#   35345 AGK-D-H3K27AC.bed|AGK-R-H3K27AC.bed
#   
  
  # get the sample ID from the overlap category
  SampleID<-venn_summary[["Name"]][grep("|",venn_summary[["Name"]],fixed = TRUE)]
  SampleID<-gsub(pattern = "^(.*)-[DR]-(.*).bed\\|.*$",replacement = "\\1-\\2",x = SampleID,perl = TRUE)
  # PROBLEM: If no overlap category present, then need to get SampleID from first entry in Names
  if(identical(SampleID, character(0))) SampleID <- gsub(pattern = "^(.*)-[DR]-(.*).bed.*$",replacement = "\\1-\\2",x = venn_summary[["Name"]][1],perl = TRUE)
  
  # need to rename the entries in the Name column
  venn_summary[["Name"]]<-sapply(X =venn_summary[["Name"]],
                                 FUN = function(z){
                                   if ( ! grepl("|",z,fixed = TRUE)) {
                                     # if its not a overlap entry change the name to D or R
                                     gsub(pattern = "^.*-([DR])-.*$",replacement = "\\1",x = z,perl = TRUE)
                                   } else if (grepl("|",z,fixed = TRUE)) {
                                     # if its an overlap entry change to 'overlap'; only overlap entries have |
                                     # gsub(pattern = "^(.*)-[DR]-(.*).bed\\|.*$",replacement = "\\1-\\2",x = z,perl = TRUE)
                                     gsub(pattern = "^.*$",replacement = "Overlap",x = z,perl = TRUE)
                                   }
                                 })
  
  
  
  # get the names to keep
  t_table_names<-venn_summary[["Name"]]
  
  # transpose the table
  t_table<-as.data.frame(t(venn_summary[["Total"]]))
  
  # set the colnames and rownames again
  colnames(t_table)<-t_table_names
  rownames(t_table)<-SampleID
  
  # check if a column is missing; if so, create it with value 0
  for(column in c("Overlap","D","R")){
    if( ! column %in% colnames(t_table)) t_table <- setNames(cbind(t_table,0), c(names(t_table), column))
  }
  
  # reverse the order of the columns so it comes out D Overlap R
  t_table <- t_table[,sort(names(t_table),decreasing = FALSE)]
  
  return(t_table)
  # example:
#   D Overlap     R
#   ZNK-H3K9ME3     67    1213 84793
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 





# ~~~~~ RUN TABLE FUNCTION ON ARGS ~~~~~~~ #
# sort the arg; this is sorting by the alphabetical on the filepaths passed
args<-args[order(args,decreasing = TRUE)]

# apply the function to all the venn.txt entries and create a single data frame out of them
overlap_df<-do.call(rbind, lapply(args, FUN = get_venn_entry))
cat("Peak overlap data frame is:\n"); overlap_df
# example:
# Peak overlap data frame is:
#   D Overlap     R
# ZNK-H3K9ME3     67    1213 84793
# ZNK-H3K9AC     651   10310 10174
# ZNK-H3K4ME3    468   11133  4090



# caluclate the percentage of the total for each entry
overlap_pcnt<-prop.table(as.matrix(overlap_df),1)
cat("Peak overlap percentage data frame is:\n"); overlap_pcnt
# example:
# Peak overlap percentage data frame is:
#   D      Overlap           R
# ZNK-H3K9ME3  0.0007784090 1.409269e-02 0.985128902
# ZNK-H3K9AC   0.0308019872 4.878164e-01 0.481381595
# ZNK-H3K4ME3  0.0298260149 7.095150e-01 0.260658976
 

# save the tables to the current dir
write.table(x = overlap_df,file = paste0(plot_filename,"_all_peak_overlap.tsv"),quote = FALSE,sep = "\t",row.names = rownames(overlap_df))
write.table(x = overlap_pcnt,file = paste0(plot_filename,"_all_peak_overlap_pcnt.tsv"),quote = FALSE,sep = "\t",row.names = rownames(overlap_pcnt))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 





# ~~~~~ SET UP THE PLOTS ~~~~~~~ #
# plot colors
BARPLOT_COLORS<-c("red","blue","green")

# setup the matrix for the plot layout
plot_layout_matrix<-structure(c(1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 1L, 2L, 
                                2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 
                                3L, 3L, 3L, 3L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L),
                              .Dim = c(11L,4L), 
                              .Dimnames = list(NULL, c("V1", "V2", "V3", "V4")))

# start the PDF output
pdf(file = paste0(plot_filename,"_peak_overlap_barplots.pdf"),width = 8,height = 10)
# setup the panel layout
layout(plot_layout_matrix)
# need to set this for some reason
par(mar=c(1,0,3,0))
# call blank plot to fill the first panel
plot(1,type='n',axes=FALSE,xlab="",ylab="",main = "Peak Overlap",cex.main=1) 
# set up the Legend in the first panel
legend("bottom",legend=colnames(overlap_df),fill=BARPLOT_COLORS,bty = "n",ncol=length(BARPLOT_COLORS),cex=1.0)
# plot margins # c(bottom, left, top, right) # default is c(5, 4, 4, 2) + 0.1
par(mar=c(5,8,0,3)+ 0.1) 

# create barplot for the two tables
barplot(t(overlap_df),cex.names = 0.7,
        horiz = T,
        col=BARPLOT_COLORS,
        border=NA,
        las=1,
        # cex.names=Names_scale,
        xlab="Number of peaks",
        space=0.6
) 

barplot(t(overlap_pcnt),cex.names = 0.7,
        horiz = T,
        col=BARPLOT_COLORS,
        border=NA,
        las=1,
        # cex.names=Names_scale,
        xlab="Percentage of peaks",
        space=0.6
) 
dev.off()


