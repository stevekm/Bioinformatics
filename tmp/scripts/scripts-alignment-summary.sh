#!/bin/bash
# USAGE: alignment_summary_stats_plots.sh /path/to/project/chipseq-analysis /path/to/project/chipseq-analysis_output
# this script will parse stats sheet to get numbers of reads
# for all samples in the project
# write all the values and sample names into a CSV
# then call an R script to do further processing of the CSV table, and create barplots to visualize the results

# $0 path to script
# $1 /path/to/project/chipseq-analysis
# $2 /path/to/output_dir 
BARPLOT_SCRIPT=~/pipeline-master/code/scripts-alignment-summary.R

##
##  Get the alignment stats, save to a table
##
##
# sample sheet is here:
FASTQ_SAMPLE_SHEET=${1}/alignments/inputs/sample-sheet.tsv
# get the sample names
SAMPLE_BASENAMES=$(cat $FASTQ_SAMPLE_SHEET | cut -d $'\t' -f2 | tail -n +2 )

# bam file dirs here
BAM_DIR=${1}/alignments/results

# the output dir for the alignment summary stats sheet
# TABLE_OUTPUT_DIR=${1}/report/report_development/alignment_summary
TABLE_OUTPUT_DIR=${2}/report_development/alignment_summary

# make the dir
mkdir -p $TABLE_OUTPUT_DIR
# the output table file
SUMMARY_CSV=${TABLE_OUTPUT_DIR}/summary_alignment.csv

# create a header for the table; use a tmp object
printf "Sample,Total Reads,Aligned Reads,De-duplicated alignments\n" > ${SUMMARY_CSV}.tmp
# cat ${SUMMARY_CSV}.tmp

for i in $(find $BAM_DIR -type f -name "stats.tsv"); do
#   echo $i; 
#   cat $i;
#   echo ""
  TOTAL_READS=$(grep "Total reads" $i | cut -d $'\t' -f2);
#   echo "  Total reads is $TOTAL_READS"
#   
  ALGN_READS=$(grep "Aligned reads" $i | cut -d $'\t' -f2);
#   echo "  Aligned reads is $ALGN_READS"
#   
  DDUP_ALGN=$(grep "De-duplicated alignments" $i | cut -d $'\t' -f2)
#   echo "  Dedup aligns is $DDUP_ALGN"
  SAMPLE_ID=$(dirname $i)
  SAMPLE_ID=$(basename $SAMPLE_ID)
  # echo $TOTAL_READS;
  printf "${SAMPLE_ID},${TOTAL_READS},${ALGN_READS},${DDUP_ALGN}\n" >> ${SUMMARY_CSV}.tmp;
  # echo ""
  done

# for i in $SAMPLE_BASENAMES; do 
#   TOTAL_READS=$(grep "Total reads" $BAM_DIR/align.$i/stats.tsv | cut -d $'\t' -f2);
#   ALGN_READS=$(grep "Aligned reads" $BAM_DIR/align.$i/stats.tsv | cut -d $'\t' -f2);
#   DDUP_ALGN=$(grep "De-duplicated alignments" $BAM_DIR/align.$i/stats.tsv | cut -d $'\t' -f2)
#   # echo $TOTAL_READS;
#   printf "${i},${TOTAL_READS},${ALGN_READS},${DDUP_ALGN}\n" >> ${SUMMARY_CSV}.tmp;
#   done

# sort combined file ignoring header and remove extra header lines
(head -n 1 ${SUMMARY_CSV}.tmp) > $SUMMARY_CSV
(tail -n +2 ${SUMMARY_CSV}.tmp | sort --unique ) >> $SUMMARY_CSV

# cat $SUMMARY_CSV
rm -f ${SUMMARY_CSV}.tmp


##
## Pass the alignement stats table to the R plot script
##
##
# 
$BARPLOT_SCRIPT ${SUMMARY_CSV} ${TABLE_OUTPUT_DIR}
