#!/bin/bash
set -x

## USAGE: peak_overlap_HOMER_mergePeak_pipeline.sh
## this pipeline script will parse our ChIP-Seq analysis results
## and run HOMER mergePeaks for each set of samples,
## then create Venn Diagrams & UpSet plots (individually), and bar plots (aggregate) to visualize the results


# ~~~~ COMMON ITEMS AND PLACES~~~~~ #
# common items 
module load homer/v4.6
# out Gencode TSS regions bed file
gen_bed="/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2016-03-31/project_data2/gencode.v19.annotation_TSSd500_10kbp.bed"
# dir with the peak calling results
pipeline_peaks_dir="/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2016-03-31/pipeline/peaks"
# outdir for all overlap sets
main_outdir="/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2016-03-31/project_notes/peak_overlap/manuscript_overlaps_report/peak_overlaps"
mkdir -p "$main_outdir"
cd "$main_outdir"

# ~~~~~~~~~~~~~~~~~~~~~~ #



# ~~~~ CUSTOM FUNCTIONS ~~~~~ #
# custom functions for the script
make_multi_venn()
{
  ## USAGE: <sampleID> venn.txt
  ## this function will make the plot for the HOMER venn file
  func_script="/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2016-03-31/project_notes/code/multi_peaks_Venn.R"
  func_ID="$1"
  func_file="$2"
  chmod +x "$func_script"
  echo "Making venn diagram.."
  module switch r/3.2.0
  "$func_script" "$func_ID" "$func_file"
}

make_UpSetR()
{
  ## USAGE: <sampleID> venn.txt
  ## this function will make the plot for the HOMER venn file
  func_script="/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2016-03-31/project_notes/code/multi_peaks_UpSet_plot.R"
  func_ID="$1"
  func_file="$2"
  chmod +x "$func_script"
  echo "Making UpSetR plot.."
  module switch r/3.3.0
  "$func_script" "$func_ID" "$func_file"
}
# ~~~~~~~~~~~~~~~~~~~~~~ #


# ~~~~ MERGEPEAKS PIPELINE ~~~~~ #
# make array items to hold the outdirs, peak files sets
unset outdir_array; unset files_array
declare -A outdir_array
declare -A files_array
declare -A marks_array
declare -a index_array
# was gonna do more stuff with these arrays but ended up not using them.. I think

# get the indexes from the branch names for the peaks
index_array=($(find ${pipeline_peaks_dir}/results -maxdepth 1 -mindepth 1 -type d ! -name ".db" -exec basename {} \; | sort | tr '\n' ' '))

# test the array..
for key in ${!index_array[@]}; do
( # sooo sooo much parallells ! 
  echo "key is ${key}"
  # get the key value
  key_value=${index_array[${key}]}
  
  # get the information for the files to be processed
  # # list of files
  # files_array["${key_value}"]=$(find ${pipeline_peaks_dir}/results/ -type f -name "peaks.bed" -path "*/${key_value}/*")
  # tmp_files=$(find ${pipeline_peaks_dir}/results/ -type f -name "peaks.bed" -path "*/${key_value}/*" | tr '\n' ' ')
  
  # # list of prefixes to make into outdir subdirs
  # outdir_array["${key_value}"]=$(find ${pipeline_peaks_dir}/results -type f -name "peaks.bed" -path "*/${key_value}/*" -exec dirname {} \; | sed -e 's|^\(.*/\)\([[:upper:]].*\)$|\2|' -e 's|^\([[:upper:]]..\).*$|\1|' | sort -u | tr '\n' ' ')
  tmp_outdir_parent=$(find ${pipeline_peaks_dir}/results -type f -name "peaks.bed" -path "*/${key_value}/*" -exec dirname {} \; | sed -e 's|^\(.*/\)\([[:upper:]].*\)$|\2|' -e 's|^\([[:upper:]]..\).*$|\1|' | sort -u | tr '\n' ' ')
  
  # list of the histone marks alone
  # marks_array["${key_value}"]=$(find ${pipeline_peaks_dir}/results -type f -name "peaks.bed" -exec dirname {} \; | sed -e 's|^\(.*/\)\([[:upper:]].*\)$|\2|' -e 's|^\([[:upper:]]..\)-[DR]-\(.*\)$|\2|' | sort -u | tr '\n' ' ')
  tmp_marks=$(find ${pipeline_peaks_dir}/results -type f -name "peaks.bed" -exec dirname {} \; | sed -e 's|^\(.*/\)\([[:upper:]].*\)$|\2|' -e 's|^\([[:upper:]]..\)-[DR]-\(.*\)$|\2|' | sort -u | tr '\n' ' ')
  
  
  # iterate over the outdirs
  for sample in $tmp_outdir_parent; do
    # iterate over the marks
    for mark  in $tmp_marks; do
      # make the outdir path
      outdir_path="${main_outdir}/${sample}/${mark}/${key_value}"
      mkdir -p "$outdir_path"
      tmp_logdir="${outdir_path}/logs"; mkdir -p "$tmp_logdir"
      
      # get the corresponding files; D and R both
      tmp_files=$(find ${pipeline_peaks_dir}/results/ -type f -name "peaks.bed" -path "*/${key_value}/*" -path "*${sample}*" -path "*${mark}*")
      num_files=$(echo "$tmp_files" | wc -l)
      tmp_files=$(echo "$tmp_files" | tr '\n' ' ')
      # tmp_file_IDs=$(find ${pipeline_peaks_dir}/results/ -type f -name "peaks.bed" -path "*/${key_value}/*" -path "*${sample}*" -path "*${mark}*" -exec dirname {} \; | sed -e 's|^\(.*/\)\([[:upper:]].*\)$|\2|')
      
      # if two files were found, do the D vs R overlap
      # # make the symlinks to the peaks files
      if [ "$num_files" == "2" ]; then
        for i in $tmp_files; 
          do echo "$i"; 
          tmp_file_ID="$(basename $(dirname $i))"; 
          echo "$tmp_file_ID"; 
          ln -fs "$i" "${outdir_path}/${tmp_file_ID}.bed"; 
        done
      fi
      # # run the mergePeaks and make plots
#       [ "$num_files" == "2" ] && qsub -o :${tmp_logdir}/ -e :${tmp_logdir}/ -N "${sample}-${mark}" <<E0F
#         module load homer/v4.6
#         cd $outdir_path
#         mergePeaks ${sample}-[DR]-${mark}.bed -prefix mergepeaks -venn venn.txt -matrix matrix.txt
#         make_multi_venn "${sample}-${mark}" venn.txt
#         make_UpSetR "${sample}-${mark}" venn.txt
#       
# E0F

      if [ "$num_files" == "2" ]; then
        (
          cd $outdir_path; 
          # if not already done, run mergePeaks and make the plots
          [ ! -f venn.txt ] && mergePeaks ${sample}-[DR]-${mark}.bed -prefix mergepeaks -venn venn.txt -matrix matrix.txt &&  make_multi_venn "${sample}-${mark}" venn.txt && make_UpSetR "${sample}-${mark}" venn.txt
          
        ) &
      
    done
    
    
    
  done
  
  
  
) &   
done

# THE END
exit
# ~~~~~~~~~~~~~~~~~~~~~~ #


# ~~~~ POST-PROCESSING ~~~~~ #
# aggregate all of the venn.txt for each analysis branch, run the R aggregagtor & barplot script on them
# run this manually when the above is all finished ; have to wait for all jobs to complete...
tmp_script="/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2016-03-31/project_notes/code/peaks_overlap_aggregatoR.R"
main_outdir="/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2016-03-31/project_notes/peak_overlap/manuscript_overlaps_report/peak_overlaps"
cd $main_outdir
# dir with the peak calling results
pipeline_peaks_dir="/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2016-03-31/pipeline/peaks"

declare -a index_array

# get the indexes from the branch names for the peaks
index_array=($(find ${pipeline_peaks_dir}/results -maxdepth 1 -mindepth 1 -type d ! -name ".db" -exec basename {} \; | sort | tr '\n' ' '))

# make the aggregate barplots off the venn.txt
for key in ${!index_array[@]}; do
  echo "key is ${key}"

  # get the key value
  key_value=${index_array[${key}]}

  # find all the venn.txt files for the branch, pass to the peaks aggregator & barplot  
  find "$main_outdir" -path "*/${key_value}/*" -name "venn.txt" -print0 | xargs -0 $tmp_script
  # NOTE: this will make the plots & table outputs in the PWD !! 
done
# ~~~~~~~~~~~~~~~~~~~~~~ #




