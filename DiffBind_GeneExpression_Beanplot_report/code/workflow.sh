#!/bin/bash
set -x

# this pipeline will:
# set up the directory parsing for the diffbind branches
# start with one branch only.. 
# run the diffbind vs gene expression script for all of them

# path to this script
# /ifs/home/kellys04/projects/SmithLab_ChIpSeq_2016-12-31/project_notes/diffferential_gene_expression/diff_expr_report/report/workflow.sh


# path to all the diffbind results
# just start with one analysis params branch for now, there are several
diffbind_results_dir="/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2016-12-31/pipeline/diffbind/results_original"
# path to the output directory
main_outdir="/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2016-12-31/project_notes/diffferential_gene_expression/diff_expr_report/analysis_output_beanplots_all_samples2"

# path to gene expression file to use for ALL comparisons..
gene_expr_file="/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2016-12-31/project_notes/diffferential_gene_expression/gene_expression.tsv"

# path to the script that compares the gene expr to diffbind
diff_expr_script="/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2016-12-31/project_notes/diffferential_gene_expression/diff_expr_report/code/diffbind_geneExpression_analysis_3.R"
## USAGE: diffbind_geneExpression_analysis.R /path/to/outdir /path/to/microarray_gene_expression.tsv /path/to/diffbind.csv <HistMark> <analysis_params_branch>
## OUTPUT: <outdir>/<SampleID>/<HistMark>/<analysis_params_branch>/[paired_boxplot.pdf genes_up.tsv genes_down.tsv]
## e.g. <outdir>/ZNK/INPUT/peaks.by_sample.macs_broad-cutoff0025
chmod +x $diff_expr_script


# get the path to all the diffbind files to use.. 
diff_files=$(find $diffbind_results_dir -type f -name "diff_bind.*" -name "*.csv")

# get the items for the script from each of the files found
for i in $diff_files; do
	# echo "$i"

	# the diffbind file found
	diff_csv_path="$i"
	echo "diff_csv_path is $diff_csv_path"

	# get the histone mark from the path
	tmp_mark="${diff_csv_path#*/align.by_sample.bowtie2/}"
	tmp_mark="${tmp_mark%/*}"
	echo "tmp_mark is $tmp_mark"

	# get the peak calling branch
	tmp_branch="${diff_csv_path#*/diffbind.by_chip.status/}"
	tmp_branch="${tmp_branch%%/*}"
	echo "tmp_branch is $tmp_branch"

	# get the file basename
	tmp_basename="$(basename $diff_csv_path)"
	tmp_basename="${tmp_basename//.csv/}"
	echo "tmp_basename is $tmp_basename"

	# concatenate the filename, and the branch
	tmp_cat_branch_basename="${tmp_branch}_${tmp_basename}"
	echo "tmp_cat_branch_basename is $tmp_cat_branch_basename"
	
	# make the output directory
	tmp_outdir="${main_outdir}/${tmp_cat_branch_basename}/${tmp_mark}"
	mkdir -p "$tmp_outdir"

	# echo "Command is:"
	# echo "Rscript --slave --no-save --no-restore $diff_expr_script $main_outdir $gene_expr_file $diff_csv_path $tmp_mark $tmp_cat_branch_basename"

	# pass the args to the script
	# module switch r/3.3.0
	# (
	# Rscript --slave --no-save --no-restore $diff_expr_script "$main_outdir" "$gene_expr_file" "$diff_csv_path" "$tmp_mark" "$tmp_cat_branch_basename"
	# ) & # run these all at once; be careful with this!
	# ) 1>&2 # write the stdout to stderr ; use this for logging later

	# submit job to qsub instead
	tmp_logdir="${tmp_outdir}/logs"; mkdir -p "$tmp_logdir"
	tmp_outdir="$tmp_outdir"; cd "$tmp_outdir"
	module unload r
	module load r/3.3.0
	qsub -b y -wd "$tmp_outdir" -o :${tmp_logdir}/ -e :${tmp_logdir}/ -pe threaded 1 "$diff_expr_script" "$tmp_outdir" "$gene_expr_file" "$diff_csv_path" "$tmp_mark" "$tmp_cat_branch_basename"


done


