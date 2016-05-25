#!/bin/bash
# this script will read the number of lines for each BED file and save it into a giant table with sample names, etc..


#~~~ Shell Options ~~~~#
# remember whether extglob was originally set, so we know whether to unset it
shopt -q extglob; extglob_set=$?
# set extglob if it wasn't originally set.
((extglob_set)) && shopt -s extglob
# Note, 0 (true) from shopt -q is "false" in a math context.

shopt -q nullglob; nullglob_set=$?
((nullglob_set)) && shopt -s nullglob

shopt -q globstar; globstar_set=$?
((globstar_set)) && shopt -s globstar

#~~~~~~~~~~~~~~~~~~#

# get the db.tsv files

# where the analysis was done
Peaks_dir="/home/steve/projects/ChIpSeq_2016/pipeline/peaks"

# the tsv file with all the params and items for the analysis
PeakResults_db="/home/steve/projects/ChIpSeq_2016/pipeline/peaks/results/.db/db.tsv"

# where to output the stats table, etc
PeakStats_dir="/home/steve/projects/ChIpSeq_2016/project_notes/peaks_stats"
# table file to make
Peakstats_table=$PeakStats_dir/stats.tsv
# makethe file
touch $Peakstats_table


# ~~~~~~ some misc scrap codes ~~~~~~~~~~~ #
# tail -n +2 $PeakResults_db | awk -F $'\t' 'BEGIN { OFS = FS } ; { print $1,$2 }'
# tail -n +2 $PeakResults_db | awk -F $'\t' 'BEGIN { OFS = FS } ; { print $NF }'

# head -n 2 pipeline/peaks/results/.db/db.tsv
# command inp-branch  inp-obj-var inp-object  out-branch  out-dir out-obj-var out-object  out-prefix  params  script  split-value split-var status
# ./code/code.main/scripts-qsub-wrapper 1 ./code/chipseq-peaks.tcsh results/peaks.by_sample.macs_broad/align.by_sample.bowtie2/AGK-D-H3K27AC params/params.macs_broad.tcsh inpdirs/align/results/align.by_sample.bowtie2 'AGK-D-H3K27AC'  inpdirs/align/results/align.by_sample.bowtie2 sample  'AGK-D-H3K27AC' results/peaks.by_sample.macs_broad/align.by_sample.bowtie2  results/peaks.by_sample.macs_broad/align.by_sample.bowtie2/AGK-D-H3K27AC  sample  AGK-D-H3K27AC results/peaks params/params.macs_broad.tcsh ./code/code.main/scripts-qsub-wrapper 1 ./code/chipseq-peaks.tcsh     up-to-date

head $PeakResults_db | awk -F $'\t' 'BEGIN { OFS = FS } ; { print $10 }' # didn't end up using this
head -n 1 $PeakResults_db  # didn't end up using this either
# command inp-branch  inp-obj-var inp-object  out-branch  out-dir out-obj-var out-object  out-prefix  params  script  split-value split-var status
# out-branch # 5 <- type of analysis
# out-dir # 6
# out-obj-var # 7
# out-object # 8 <- sample name
# out-prefix # 9
# ~~~~~~ some misc scrap codes ~~~~~~~~~~~ #



SampleIDs=$(tail -n +2 $PeakResults_db | cut -f8 | sort -u)
SampleNames=$(tail -n +2 $PeakResults_db | cut -f8 | cut -f1 -d "-" | sort -u) # didn't need this one either
Branch=$(tail -n +2 $PeakResults_db | cut -f5 | sort -u)

# wipe the output file if it has contents
echo -ne "" > $Peakstats_table
# print the header
echo -e "Patient\tStatus\tChIP\tSample\tNumPeaks\tGrouping\tFile"> $Peakstats_table
for i in $Branch; do
  # flatten the branch name
  tmp_Branch=$(echo "$i" | tr "/" "-")

  # make output dir
  #mkdir -p $PeakStats_dir/$tmp_Branch
  
  #echo -e "\n$i\n"
  
  # get the bed files for the branch
  BEDFILES=($Peaks_dir/$i/**/peaks.bed)
  # echo "${BEDFILES[@]}" | tr " " "\n"
  # for q in $SampleNames; do
  for q in $SampleIDs; do
    tmp_SampleName=$(echo $q | cut -f1 -d "-")
    tmp_Status=$(echo $q | cut -f2 -d "-")
    tmp_ChIP=$(echo $q | cut -f3 -d "-")
    # make outdir
    #mkdir -p $PeakStats_dir/$tmp_Branch/$q
    # echo $q
    
    # find the bed files that contain the given sample name
    BEDmatches=$(echo "${BEDFILES[@]}" | tr " " "\n" | grep "$q")
    
    # count the number of lines in the bed file = peaks
    NumPeaks=$(cat $BEDmatches | wc -l )
    
    # the parent dir is the sample name full ID, get it
    # tmp_dirname=""
    
    
    # printf "$BEDmatches\t$i\t$tmp_SampleName\t$q\t$NumPeaks\t\n"
    # printf "$BEDmatches\t$i\t$tmp_SampleName\t$q\t$NumPeaks\t\n" >> $Peakstats_table

    # print the entries into the table
    printf "$tmp_SampleName\t$tmp_Status\t$tmp_ChIP\t$q\t$NumPeaks\t$i\t$BEDmatches\t\n" >> $Peakstats_table
  done
done

# looks like this:
#Patient Status  ChIP  Sample  NumPeaks  Grouping  File
#AGK D H3K27AC AGK-D-H3K27AC 14995 results/peaks.by_group.macs_broad/align.by_sample.bowtie2 /home/steve/projects/ChIpSeq_2016/pipeline/peaks/results/peaks.by_group.macs_broad/align.by_sample.bowtie2/AGK-D-H3K27AC/peaks.bed
#AGK D H3K27ME3  AGK-D-H3K27ME3  0 results/peaks.by_group.macs_broad/align.by_sample.bowtie2 /home/steve/projects/ChIpSeq_2016/pipeline/peaks/results/peaks.by_group.macs_broad/align.by_sample.bowtie2/AGK-D-H3K27ME3/peaks.bed
#AGK D H3K4ME3 AGK-D-H3K4ME3 12730 results/peaks.by_group.macs_broad/align.by_sample.bowtie2 /home/steve/projects/ChIpSeq_2016/pipeline/peaks/results/peaks.by_group.macs_broad/align.by_sample.bowtie2/AGK-D-H3K4ME3/peaks.bed





#~~~~~~~~~~~~~~~~~~#
# unset globs if it wasn't originally set
((extglob_set)) && shopt -u extglob
((nullglob_set)) && shopt -u nullglob
((globstar_set)) && shopt -u globstar
