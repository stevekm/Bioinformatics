#!/bin/bash
# find the output bams from the pipeline, make a quick sample sheet, then use that to run MACS1 on them

# locations
pipeline_sample_sheet="/ifs/home/kellys04/projects/SmithLab_ChIPSeq_RNASeq_2016_06_01/inputs/sample-sheet.tsv"
align_results_dir="/ifs/home/kellys04/projects/SmithLab_ChIPSeq_RNASeq_2016_06_01/pipeline/align/results"
analysis_dir="/ifs/home/kellys04/projects/SmithLab_ChIPSeq_RNASeq_2016_06_01/project_notes/Combined_ChIP_RNA_Seq"
analysis_outdir="${analysis_dir}/analysis_output"
align_sample_sheet="${analysis_dir}/alignment_sample-sheet.tsv"

# get the bam files from the pipeline for looping over
bam_files=$(find $align_results_dir -name "*.bam" | tr '\n' ' ')


# make an empty file; clear the file if it already exists
echo -n "" > $align_sample_sheet

# iterate over the bams
for i in $bam_files; do
  tmp_bam="$i"
  tmp_ID="$(basename $(dirname $tmp_bam))"
  echo $tmp_ID
  echo $tmp_bam
  
  # write to the file
  echo -e "${tmp_ID}\t${tmp_bam}" >> $align_sample_sheet
done

# #
# manually added third column in Excel for the Input BAM !!
# # 

# new sample sheet:
align_sample_sheet="${analysis_dir}/alignment_sample-sheet2.tsv"

# read in the sample sheet, iterate over each line
cat $align_sample_sheet | while read i; do
# make sure there is an entry on that line !
if [[ ! -z "$i" ]]; then
		tmp_sample=$(echo "$i" | cut -f1)
		tmp_sample_bam=$(echo "$i" | cut -f2)
		tmp_input_bam=$(echo "$i" | cut -f3)
		# echo "$tmp_input_bam"
		
		# make sure the input entry is not NA e.g. the Input sample
		# [[ ! $tmp_input_bam == "NA" ]]  && 
		if [[ ! $tmp_input_bam == "NA" ]]; then
			echo "tmp_sample is $tmp_sample"
			echo "tmp_sample_bam is $tmp_sample_bam"
			echo "tmp_input_bam is $tmp_input_bam"
			
			# make outdir
			sample_outdir="${analysis_outdir}/macs14/${tmp_sample}"
			tmp_logdir="${sample_outdir}/logs"
			mkdir -p "${tmp_logdir}"
			
			qsub -wd $sample_outdir -o :${tmp_logdir}/ -e :${tmp_logdir}/ -N "$tmp_sample" <<E0F1
				#!/bin/bash
				set -x
				module unload python
				module load macs/1.4.2
				echo "test"
				echo "pwd is"
				pwd
				echo "tmp_sample is $tmp_sample"
				echo "tmp_sample_bam is $tmp_sample_bam"
				echo "tmp_input_bam is $tmp_input_bam"
				macs14  --format=BAM --gsize=hs --bdg --single-profile --diag --name=$tmp_sample -t $tmp_sample_bam -c $tmp_input_bam
			
E0F1
		fi
	fi
done
