#!/bin/bash
# set -x

# this workflow will run bowtie on the input fastq files
# then run MACS1 and generate peaks and bedGraphs



ProjDir="/ifs/home/kellys04/projects/SmithLab_ChIPSeq_RNASeq_2016_06_01/project_notes/Combined_ChIP_RNA_Seq"
analysis_input="${ProjDir}/analysis_input"
fastq_dir="${analysis_input}/fastq2"
analysis_outdir="${ProjDir}/analysis_output"

bam_dir="${analysis_outdir}/bam"
bam_log_dir="${bam_dir}/logs"
mkdir -p "$bam_log_dir"

peaks_dir="${analysis_outdir}/peaks"
peaks_log_dir="${peaks_dir}/logs"
mkdir -p "$peaks_log_dir"


bigwig_dir="${analysis_outdir}/bigWig"
bigwig_log_dir="${bigwig_dir}/logs"
mkdir -p "$bigwig_log_dir"


# reference genome info
GENOME="hg19"
GenomeRefDir="/local/data/iGenomes/Homo_sapiens/UCSC/${GENOME}/Sequence/Bowtie2Index/genome"




# # make a quick sample sheet to match sample names, fastq, and inputs
# setup sample sheet
# sample_sheet="${ProjDir}/bowtie_MACS_sample-sheet.tsv"
# clear the sheet
# echo -n "" $sample_sheet
# fastq_files="$(find ${fastq_dir} -name "*.fastq.gz")"
# 
# for i in $fastq_files; do
# 	tmp_fastq="$i"
# 	echo "tmp_fastq is $tmp_fastq"
# 	
# 	# get the file name
# 	tmp_ID="$(basename $tmp_fastq)"
# 	tmp_ID="${tmp_ID%.fastq.gz}"
# 	echo "tmp_ID is $tmp_ID"
# 	
# 	# write a quick sample sheet
# 	echo -e "${tmp_ID}\t${tmp_fastq}" >> $sample_sheet
# 
# done




# NEW sample sheet:
# do quick edits in Excel to match the Inputs and set groups
sample_sheet="${ProjDir}/bowtie_MACS_sample-sheet2.tsv"
# sample sheet now has header:
# SampleID        Input   Fastq   Genome

# get the items for each from the sample
tail -n +2 $sample_sheet | while read i; do
	# echo "$i"
	if [[ ! -z "$i" ]]; then
		# echo "$i"
		tmp_ID="$(echo "$i" | cut -f1)"
		echo "$tmp_ID"
		
		tmp_input="$(echo "$i" | cut -f2)"
		echo "$tmp_input"
		
		tmp_fastq="$(echo "$i" | cut -f3)"
		echo "$tmp_fastq"
		
		tmp_out_bam_path="${bam_dir}/${tmp_ID}.bam"
		echo "$tmp_out_bam_path"
		
		# submit alignment job
		# need to hard-code the threads for the heredoc...		# THREADS=${NSLOTS:=8} 
		qsub -wd $bam_dir -o :${bam_log_dir}/ -e :${bam_log_dir}/ -pe threaded 4-12 -N "$tmp_ID" <<E0F1
			#!/bin/bash
			set -x
			
			THREADS=\${NSLOTS} 
			
			module load bowtie2/2.2.6
			module unload samtools
			module load samtools/1.2.1
			
			
			
			echo "\${THREADS} is num slots"
			echo "$tmp_ID"
			echo "$tmp_input"
			echo "$tmp_fastq"
			echo "$tmp_out_bam_path"
			
			# alignment step
			if [ ! -f ${tmp_out_bam_path} ]; then
				echo "creating alignment file"
				bowtie2 -x ${GenomeRefDir} --threads 8 --local -U ${tmp_fastq} | samtools view -Sb - > ${tmp_out_bam_path}
			fi
			
			
E0F1
		
		# set up the MACS job
		# check for NA as input 
		if [[ ! $tmp_input == "NA" ]]; then
			echo "Input isnt NA"
			echo "input is $tmp_input"
			
			# check if both bam's exist
			[ -f $tmp_out_bam_path ] && echo "Sample bam exists"
			
			# find the input bam
			tmp_input_bam="$(find $bam_dir -name "${tmp_input}.bam")"
			# make sure it exists
			[[ ! -z $tmp_input_bam ]] && [ -f $tmp_input_bam ] && echo "tmp_input_bam is $tmp_input_bam"
			
			# MACS job submission
			qsub -wd $peaks_dir -o :${peaks_log_dir}/ -e :${peaks_log_dir}/ -N "$tmp_ID" <<E0F2
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
				macs14  --format=BAM --gsize=hs --bdg --single-profile --diag --name=$tmp_ID -t $tmp_out_bam_path -c $tmp_input_bam
			
E0F2
			
		fi
		
	fi
done

