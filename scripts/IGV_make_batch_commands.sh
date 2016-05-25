#!/bin/bash

# this is a basic script to read in data from ANNOVAR annotation output, and output the chromosomal coordinates into 
# a series of IGV batch script commands
# this is for taking automated screenshots / 'snapshots' in IGV of the selected loci


# file with the found mutations
Annotation_File="/home/.../projects/genomic_reporting/annotations/discovered_variants_AAmatch_unique_Varscan.csv"
# the file looks like this:
#"Chr","Start","End","Ref","Alt","Gene.refGene","AAmatch","SampleName","Program"
#"chr1",65332620,65332620,"C","T","JAK1","G307S","N.13-619.T.15-238","Varscan"
#"chr1",27105617,27105617,"C","T","ARID1A","T1743M","N.15-239.T.15-240","Varscan"
#"chr10",43615611,43615611,"G","A","RET","R897Q","N.13-619.T.15-1136","Varscan"
#"chr10",89717672,89717672,"C","T","PTEN","R233","N.15-058.T.15-060","Varscan"
#"chr11",22646800,22646800,"G","A","FANCF","A186V","N.13-619.T.15-190","Varscan"
#"chr11",108158388,108158388,"A","G","ATM","H1352R","N.13-619.T.15-237","Varscan"
#"chr11",118373862,118373862,"G","A","KMT2A","E2416K","N.13-619.T.15-237","Varscan"
#"chr11",108206576,108206576,"G","A","ATM","R2719H","N.15-236.T.15-235","Varscan"
#"chr12",56492623,56492623,"G","A","ERBB3","E925K","N.13-619.T.15-237","Varscan"

# a file to hold just the coordinates from the annotations
CoorList_file=$(dirname $Annotation_File)/IGV_chr_coor_Varscan_AAmatch.txt


# get the list of bams and their names.. for the URL for IGV
# the bam's are all hosted on our data server, accessible with http link
BamDir="/data/.../path/to/bams"
ls -1 **/*.bam > $(dirname $Annotation_File)/bam_list.txt # get a list of the bams

# make a file to hold the URLs
URLfile=$(dirname $Annotation_File)/URL_list.txt
echo -e "" > $URLfile # make the file, also could use touch but this clears the file first

# paste together the full URL for the bam files, save it into the text file
for i in $(ls $BamDir/**/*.bam); do # you're not supposed to parse the output of ls but I'm doing it anyway.
  tmp_filename=${i##$(dirname $(dirname $i))/} # strip the file path up to the last dir from the start; bam path for URL
  # echo $tmp_filename
  tmp_dirname=$(basename $(dirname $i)) # name of the sample
  # echo $tmp_dirname
  echo $tmp_dirname","'https://dataserver.com/bams/'$tmp_filename >> $URLfile
  # use awk later to parse this back out.. for the IGV commands
done

# The URL list file looks like this now:
# 
# 13-619,https://dataserver.com/bams/13-619/13-619.recalibrated.bam
# 15-058,https://dataserver.com/bams/15-058/15-058.recalibrated.bam
# 15-060,https://dataserver.com/bams/15-060/15-060.recalibrated.bam
# 15-1136,https://dataserver.com/bams/15-1136/15-1136.recalibrated.bam
# 15-1172,https://dataserver.com/bams/15-1172/15-1172.recalibrated.bam
# 15-172,https://dataserver.com/bams/15-172/15-172.recalibrated.bam
# 15-190,https://dataserver.com/bams/15-190/15-190.recalibrated.bam
# 15-235,https://dataserver.com/bams/15-235/15-235.recalibrated.bam
# 15-236,https://dataserver.com/bams/15-236/15-236.recalibrated.bam
# 15-237,https://dataserver.com/bams/15-237/15-237.recalibrated.bam
# 15-238,https://dataserver.com/bams/15-238/15-238.recalibrated.bam
# 15-239,https://dataserver.com/bams/15-239/15-239.recalibrated.bam
# 15-240,https://dataserver.com/bams/15-240/15-240.recalibrated.bam
# 15-241,https://dataserver.com/bams/15-241/15-241.recalibrated.bam
# B15-1024,https://dataserver.com/bams/B15-1024/B15-1024.recalibrated.bam
# B15-1027,https://dataserver.com/bams/B15-1027/B15-1027.recalibrated.bam


# get the coordinates in the form IGV likes, save them to the coordinates file
tail -n +2 $Annotation_File | tr -d '"' | awk -F $',' 'BEGIN { OFS = FS } ; { print $1":"$2"-"$3 }' > $CoorList_file

# the output looks like this:
#chr1:65332620-65332620
#chr1:27105617-27105617
#chr10:43615611-43615611
#chr10:89717672-89717672
#chr11:22646800-22646800
#chr11:108158388-108158388
#chr11:118373862-118373862
#chr11:108206576-108206576
#chr12:56492623-56492623
#chr12:25398281-25398281


# get the coordinates again, for iterating in the loop
CoorList=$(cat $CoorList_file) # the coordinates
# get the URLs again
URLlist=$(cat $URLfile)

# file to put the IGV commands in
# IGVCommands_file=$(dirname $Annotation_File)/IGV_snapshot_commands_Varscan_AAmatch.txt
IGVCommands_file=$(dirname $Annotation_File)/IGV_snapshot_commands_Varscan_AAmatch-no_CollapseSort.txt
echo -e "" > $IGVCommands_file

# make the commands for IGV
# for i in $CoorList; do
for i in $URLlist; do
  echo "new" >> $IGVCommands_file
  echo "snapshotDirectory /Users/steve/IGV_snapshot/"$(echo $i | awk -F $',' 'BEGIN { OFS = FS } ; { print $1 }') >> $IGVCommands_file
  echo "load "$(echo $i | awk -F $',' 'BEGIN { OFS = FS } ; { print $2 }') >> $IGVCommands_file
  echo "genome hg19" >> $IGVCommands_file
  for q in $CoorList; do
    echo "goto "$q >> $IGVCommands_file
    # echo "sort" >> $IGVCommands_file
    # echo "collapse" >> $IGVCommands_file
    echo "snapshot" >> $IGVCommands_file
  done
done


# looks like this

#new
#snapshotDirectory /Users/steve/IGV_snapshot/13-619
#load https://dataserver.com/bams/13-619/13-619.recalibrated.bam
#genome hg19
#goto chr1:65332620-65332620
#snapshot
#goto chr1:27105617-27105617
#snapshot
#goto chr10:43615611-43615611
#snapshot






# load all bams at once # this is the better way to display them !! Looks nicer I think, but the image is really really long

IGVCommands_file2=$(dirname $Annotation_File)/IGV_snapshot_commands_Varscan_AAmatch-allBAMs.txt
echo -e "" > $IGVCommands_file2
echo "new" >> $IGVCommands_file2
echo "snapshotDirectory /Users/steve/IGV_snapshot/allBAMs" >> $IGVCommands_file2
for i in $URLlist; do
  echo "load "$(echo $i | awk -F $',' 'BEGIN { OFS = FS } ; { print $2 }') >> $IGVCommands_file2
done
echo "genome hg19" >> $IGVCommands_file2
for q in $CoorList; do
  echo "goto "$q >> $IGVCommands_file2
  echo "snapshot" >> $IGVCommands_file2
done


# looks like this; truncated
#new
#snapshotDirectory /Users/steve/IGV_snapshot/allBAMs
#load https://dataserver.com/bams/13-619/13-619.recalibrated.bam
#load https://dataserver.com/bams/15-058/15-058.recalibrated.bam
#load https://dataserver.com/bams/15-060/15-060.recalibrated.bam
#load https://dataserver.com/bams/15-1136/15-1136.recalibrated.bam
#load https://dataserver.com/bams/15-1172/15-1172.recalibrated.bam
#load https://dataserver.com/bams/15-172/15-172.recalibrated.bam
#genome hg19
#goto chr1:65332620-65332620
#snapshot
#goto chr1:27105617-27105617
#snapshot
#goto chr10:43615611-43615611
#snapshot
#goto chr10:89717672-89717672
