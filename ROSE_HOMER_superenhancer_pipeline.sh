#!/bin/bash
set -x

## USAGE: ROSE_HOMER_superenhancer_pipeline.sh
## This pipeline will process previously called peaks files, perform peak merging and overlapping, and then call super enchancers with both ROSE and HOMER
## also this sets up some custom UCSC tracks but I have since made a better script for that which uses the files produced here
# http://homer.salk.edu/homer/ngs/peaks.html
# http://younglab.wi.mit.edu/super_enhancer_code.html

# ~~~ SOFTWARES ~~~~ # 
module load samtools/1.3
module load homer/v4.6
rose_dir="/ifs/home/kellys04/software/rose"
bed_to_gff_conv="/ifs/home/kellys04/software/bed_to_gff_converter.py"
# ~~~~~~~~~~~~~~~~ # 



# ~~~~ COMMON ITEMS & PLACES ~~~~ #
# out Gencode TSS regions bed filea
gen_bed="/ifs/home/kellys04/projects/SmithLab_ChIPSeq_2016-03-31/project_data2/gencode.v19.annotation_TSSd500_10kbp.bed"
# dir with the peak calling results
pipeline_peaks_dir="/ifs/home/kellys04/projects/SmithLab_ChIPSeq_2016-03-31/pipeline/peaks"
pipeline_align_dir="/ifs/home/kellys04/projects/SmithLab_ChIPSeq_2016-03-31/pipeline/align"
pipeline_homer_tagdir="/ifs/home/kellys04/projects/SmithLab_ChIPSeq_2016-03-31/pipeline/homer_tagdir"
venn_script="/ifs/home/kellys04/projects/SmithLab_ChIPSeq_2016-03-31/project_notes/code/multi_peaks_Venn.R"
UpSet_script="/ifs/home/kellys04/projects/SmithLab_ChIPSeq_2016-03-31/project_notes/code/multi_peaks_UpSet_plot.R"
# outdir for all overlap sets
main_outdir="/ifs/home/kellys04/projects/SmithLab_ChIPSeq_2016-03-31/project_notes/super_enhancer_calling/ROSE_superenhancers_2"
# a dir to hold copies of all of the peaks files produced on the external for use in UCSC
external_results_dir="/ifs/data/sequence/results/smithlab/ChIPSeq/2016-01-04/ChIP-Seq"
external_super_dir="${external_results_dir}/superenhancer"
mkdir -p "$external_super_dir"

external_bigwig_dir="/ifs/data/sequence/results/smithlab/ChIPSeq/2016-01-04/ChIP-Seq/alignment/bigwigs"
UCSC_bigwig_view='visibility=full autoScale=off alwaysZero=on maxHeightPixels=50 graphType=bar viewLimits=0:0.3'
myurlbase="https://genome.med.nyu.edu/results/"
hg19_chrom_sizes="/ifs/home/kellys04/software/hg19.chrom.sizes"
bed_to_bigbed_path="/ifs/home/kellys04/software/UCSC/bedToBigBed"
bedsort_path="/ifs/home/kellys04/software/UCSC/bedSort"
# http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/

# clear the UCSC file
echo -n "" > "${main_outdir}/UCSC_custom_track.txt.tmp"

mkdir -p "$main_outdir"
echo "main_outdir is $main_outdir"
cd "$main_outdir"
# ~~~~~~~~~~~~~~~~ # 




# ~~~ GET THE SAMPLES AND THEIR FILES ~~~~~ # 
# need the following for each sample; D and R run separately!
# # sample H3K27AC peaks.bed file
# # sample H3K27AC bam and bai files
# # sample H3K4ME3 peaks.bed file
# # gencode TSS bed file previously made
# # ref genome type e.g. hg19


# make array items to hold the outdirs, peak files sets
unset outdir_array; unset sample_peaks_files_array
declare -A outdir_array
declare -A sample_peaks_files_array

# get the peaks calling branches from the peaks results dir
peaks_branches=$( find $pipeline_peaks_dir/results/ -maxdepth 1 -type d ! -name ".db" ! -name "results" -exec basename {} \; )

echo -e "peaks_branches are:\n${peaks_branches}"

for branch in $( echo "$peaks_branches" | tr '\n' ' ' ); do 
  echo "branch is $branch"
  
  # outdir for the branch
  outdir_array[$branch]="${main_outdir}/${branch}"
  echo "outdir is ${outdir_array[$branch]}"
  mkdir -p "${outdir_array[$branch]}"
  
  # make a common dir for the peaks output per branch
  # ended up not needing this... 
  
  common_branch_peaks="${outdir_array[$branch]}_peaks_for_browser"
  mkdir -p "$common_branch_peaks"
  echo "common peak outdir is $common_branch_peaks"
  
  
  # peak files for the branch H3K27AC sample
  sample_peaks_files_array[$branch]=$(find $pipeline_peaks_dir/results/ -name "peaks.bed" -path "*/${branch}/*" -path "*/*H3K27AC*/*")
  # echo "files are ${sample_peaks_files_array[$branch]}"
  
  
  # process each of the branch H3K27AC sample files
  for sample_H3K27AC_peaks in $(echo ${sample_peaks_files_array[$branch]} | tr '\n' ' '); do
  # (
    echo -e "sample_H3K27AC_peaks is ${sample_H3K27AC_peaks}"
    
    # get the sample ID
    tmp_SampleID=$(basename $(dirname $sample_H3K27AC_peaks))
    echo "tmp_SampleID is $tmp_SampleID"
    
    # make a subdir for the output
    tmp_sample_outdir="${outdir_array[$branch]}/${tmp_SampleID}"
    mkdir -p "$tmp_sample_outdir"; cd "$tmp_sample_outdir"
    echo "tmp_sample_outdir is $tmp_sample_outdir"
    
    # make a subsubdir to contain the peaks files for merging since its gonna get messy!
    tmp_sample_peak_outdir="${tmp_sample_outdir}/peaks"
    mkdir -p "$tmp_sample_peak_outdir"
    echo "tmp_sample_peak_outdir is $tmp_sample_peak_outdir"
    
    # make logdir for qsub
    tmp_logdir="${tmp_sample_outdir}/logs"
    echo "tmp_logdir is $tmp_logdir"
    mkdir -p "$tmp_logdir"
    
    # find the sample H3K4ME3 peaks file to overlap against
    # # strip the H3K27AC from the SampleID to use for search; truncated sample ID
    tmp_SampleID_trunc=${tmp_SampleID%%-H*}
    echo "tmp_SampleID_trunc is $tmp_SampleID_trunc"
    # get the tmp_sample_H3K4ME3_peaks peaks file
    tmp_sample_H3K4ME3_peaks=$(find $pipeline_peaks_dir/results/ -name "peaks.bed" -path "*/${branch}/*" -path "*/${tmp_SampleID_trunc}-H3K4ME3/*")
    echo "tmp_sample_H3K4ME3_peaks is $tmp_sample_H3K4ME3_peaks"
    # some of the files might not exist, check the length of the find output; stop this iteration of the loop if its 0 length
    [ -z "$tmp_sample_H3K4ME3_peaks" ] && echo "No H3K4ME3 peaks found for $tmp_SampleID_trunc" && continue
    echo "tmp_sample_H3K4ME3_peaks is $tmp_sample_H3K4ME3_peaks"
    
    # get the sample name for the H3K4ME3 peaks file
    tmp_SampleID_H3K4ME3="$(basename $(dirname $tmp_sample_H3K4ME3_peaks))"
    echo "tmp_SampleID_H3K4ME3 is $tmp_SampleID_H3K4ME3"
    
    # make symlinks to these bed files in the peaks subsubdir
    ln -fs ${sample_H3K27AC_peaks} ${tmp_sample_peak_outdir}/${tmp_SampleID}.bed
    ln -fs ${tmp_sample_H3K4ME3_peaks} ${tmp_sample_peak_outdir}/${tmp_SampleID_H3K4ME3}.bed
    
    # copy over the gencode bed
    cat "$gen_bed" | cut -f1-3 > "${tmp_sample_peak_outdir}/Gencode.bed"
    
    
    
    # symlink to the H3K27AC sample bam and bai files..
    echo "Getting the alignment files and / or tagdirs"
    echo "pwd is:"
    pwd
    echo ""
    
    cd "${tmp_sample_outdir}"
    echo "check the pwd again:"
    pwd
    echo ""
    
    # rm -f "${tmp_sample_outdir}/${tmp_SampleID}.bam"
    sample_H3K27AC_bam=$(find $pipeline_align_dir/ -type f -name "alignments.bam" -path "*/*${tmp_SampleID}*/*")
    echo "sample_H3K27AC_bam is $sample_H3K27AC_bam"
    [ -z ${sample_H3K27AC_bam} ] && echo "No sample_H3K27AC_bam found; skipping" && continue # exit # && 
    [ ! -z ${sample_H3K27AC_bam} ] && ln -fs "${sample_H3K27AC_bam}" "${tmp_SampleID}.bam"
     # [ ! -f ${tmp_sample_outdir}/${tmp_SampleID}.bam ] && cp "${sample_H3K27AC_bam}" "${tmp_sample_outdir}/${tmp_SampleID}.bam"
    
    # rm -f "${tmp_sample_outdir}/${tmp_SampleID}.bam.bai"
    sample_H3K27AC_bai=$(find $pipeline_align_dir/ -type f -name "alignments.bam.bai" -path "*/*${tmp_SampleID}*/*")
    echo "sample_H3K27AC_bai is $sample_H3K27AC_bai"
    [ -z ${sample_H3K27AC_bai} ] && echo "No sample_H3K27AC_bai found; skipping" && continue #exit 
    [ ! -z ${sample_H3K27AC_bai} ] && ln -fs "${sample_H3K27AC_bai}" "${tmp_SampleID}.bam.bai"
    # [ ! -f ${tmp_sample_outdir}/${tmp_SampleID}.bam.bai ] && cp "${sample_H3K27AC_bai}" "${tmp_sample_outdir}/${tmp_SampleID}.bam.bai"
    
    # also to the INPUT alignment for this sample..
    # rm -f "${tmp_sample_outdir}/${tmp_SampleID_trunc}-INPUT.bam"
    sample_INPUT_bam=$(find $pipeline_align_dir/ -type f -name "alignments.bam" -path "*/*${tmp_SampleID_trunc}-INPUT*/*")
    echo "sample_INPUT_bam is $sample_INPUT_bam"
    [ -z ${sample_INPUT_bam} ] && echo "No sample_INPUT_bam found; skipping" && continue # exit # 
    [ ! -z ${sample_INPUT_bam} ] && ln -fs "${sample_INPUT_bam}" "${tmp_SampleID_trunc}-INPUT.bam"
    # [ ! -f ${tmp_sample_outdir}/${tmp_SampleID_trunc}-INPUT.bam ] && cp "${sample_INPUT_bam}" "${tmp_sample_outdir}/${tmp_SampleID_trunc}-INPUT.bam"
    
    # rm -f "${tmp_sample_outdir}/${tmp_SampleID}.bam.bai"
    sample_INPUT_bai=$(find $pipeline_align_dir/ -type f -name "alignments.bam.bai" -path "*/*${tmp_SampleID_trunc}-INPUT*/*")
    echo "sample_INPUT_bai is $sample_INPUT_bai"
    [ -z ${sample_INPUT_bai} ] && echo "No sample_INPUT_bai found; skipping" && continue #exit 
    [ ! -z ${sample_INPUT_bai} ] && ln -fs "${sample_INPUT_bai}" "${tmp_SampleID_trunc}-INPUT.bam.bai"
    # [ ! -f ${tmp_sample_outdir}/${tmp_SampleID}.bam.bai ] && cp "${sample_INPUT_bai}" "${tmp_sample_outdir}/${tmp_SampleID_trunc}-INPUT.bam.bai"
    
    
    # symlink to the sample and INPUT HOMER tagdirs
    # pipeline_homer_tagdir
    sample_H3K27AC_tagdir=$(find $pipeline_homer_tagdir/results/ -type d -name "*${tmp_SampleID}*")
    echo "sample_H3K27AC_tagdir is $sample_H3K27AC_tagdir"
    [ -z ${sample_H3K27AC_tagdir} ] && echo "No sample tagdir dir found; skipping" && continue
    [ ! -z ${sample_H3K27AC_tagdir} ] && ln -fs ${sample_H3K27AC_tagdir} "${tmp_SampleID}_tagdir"
    
    sample_INPUT_tagdir=$(find $pipeline_homer_tagdir/results/ -type d -name "*${tmp_SampleID_trunc}-INPUT*")
    echo "sample_INPUT_tagdir is $sample_INPUT_tagdir"
    [ -z ${sample_INPUT_tagdir} ] && echo "No input tagdir found; skipping" && continue
    [ ! -z ${sample_INPUT_tagdir} ] && ln -fs ${sample_INPUT_tagdir} "${tmp_SampleID_trunc}-INPUT_tagdir"
    
    # symlink to all the files in the rose_dir because thats how ROSE rolls...
    # this is the directory that contains all the program files for ROSE as downloaded from the Young lab website
    ln -fs ${rose_dir}/* ${tmp_sample_outdir}/
    
    # need to get the matching bigWig from the external dir..
    # /ifs/data/sequence/results/smithlab/ChIPSeq/2016-01-04/ChIP-Seq/alignment/bigwigs
    sample_bigWig=$(find "$external_bigwig_dir" -type f -name "*${tmp_SampleID}*")
    # [ -z ]
    
    # run HOMER mergePeaks; cd to that dir for it though..
    # tmp_pwd=$PWD; 
    
    # this can be submitted to qsub from here on as a HereDoc:
    # need to run this part of the script on the HPC cluster! Run it all in parallel its very resouce intensive!
    # qsub -o :${tmp_logdir}/ -e :${tmp_logdir}/ -pe threaded 1 -l h_vmem=15G -l mem_free=15G -l mem_token=15G <<E0F
    qsub -o :${tmp_logdir}/ -e :${tmp_logdir}/ -N ${tmp_SampleID}_${branch} <<E0F
      #!/bin/bash # uncomment this for qsub .. or maybe not
      set -x
      module load samtools/1.3
      module load homer/v4.6
      
      cd "${tmp_sample_peak_outdir}"
      # do the mergePeaks & plot the results
      [ ! -f venn.txt ] && mergePeaks "${tmp_SampleID}.bed" "${tmp_SampleID_H3K4ME3}.bed" Gencode.bed -prefix mergepeaks -venn venn.txt -matrix matrix.txt
      [ ! -f ${tmp_SampleID}_venn.pdf ] && module switch r r/3.2.0 && Rscript /ifs/home/kellys04/projects/myBioinformatics/HOMER_mergePeaks_multiVenn/multi_peaks_Venn.R "${tmp_SampleID}" venn.txt
      [ ! -f ${tmp_SampleID}_UpSetR_plot.pdf ] && module switch r r/3.3.0 && Rscript /ifs/home/kellys04/projects/myBioinformatics/HOMER_mergePeaks_venn_UpSetR/multi_peaks_UpSet_plot.R "${tmp_SampleID}" venn.txt
      # HOMER gives messed up BED outputs have to clean up HOMER's mess and reformat the bed files..
      # path to the peaks unique to the H3K27AC
      # tmp_ROSE_bed="$( readlink -m mergepeaks_${tmp_SampleID}.bed )"; echo "tmp_ROSE_bed is $tmp_ROSE_bed"; wc -l $tmp_ROSE_bed; echo ""
      # tmp_ROSE_bed="$( readlink -m ${tmp_sample_peak_outdir}/mergepeaks_${tmp_SampleID}.bed )"; echo "tmp_ROSE_bed is $tmp_ROSE_bed"; wc -l $tmp_ROSE_bed; echo ""
      
      # reformat back to vanilla bed file format..
      tail -n +2 mergepeaks_${tmp_SampleID}.bed | cut -f2,3,4,9 > mergepeaks_${tmp_SampleID}.bed_goodbed
      
      # convert the bed file to gff; remove the dumb header
      python ${bed_to_gff_conv} "mergepeaks_${tmp_SampleID}.bed_goodbed" "mergepeaks_${tmp_SampleID}.bed_goodbed.gff" ; sed -i 1,3d "mergepeaks_${tmp_SampleID}.bed_goodbed.gff"
      # clean some stupid suff the program left in the gff
      # I broke this step last time and ROSE broke and had massive memory leak watch out don't mess this part up !! note the inserted tab in the sed command very important
      cat "mergepeaks_${tmp_SampleID}.bed_goodbed.gff" | sed -e 's/bed2gff[[:space:]]/\t/g' > "mergepeaks_${tmp_SampleID}.bed_goodbed.gff.better.gff"
      #  sed -e 's/bed2gff[[:space:]]/\t/g'
      cd "$tmp_sample_outdir"
      
      
      
      # run the ROSE 
      echo -e "\npwd is $PWD \n"
      pwd
      ls -l *
      
      echo -e "The command for ROSE is:\n"
      echo "python ./ROSE_main.py -g HG19 -i ${tmp_sample_peak_outdir}/mergepeaks_${tmp_SampleID}.bed_goodbed.gff.better.gff -r ${tmp_SampleID}.bam -c ${tmp_SampleID_trunc}-INPUT.bam -o ROSE_super_enhancer -t 2500"
      echo "calling ROSE:"
      python ./ROSE_main.py
      echo "Writing the command for ROSE to the file.."
      echo -e '#!/bin/bash' > ROSE_command.sh
      echo "set -x; cd $tmp_sample_outdir; module load samtools; python ./ROSE_main.py -g HG19 -i ${tmp_sample_peak_outdir}/mergepeaks_${tmp_SampleID}.bed_goodbed.gff.better.gff -r ${tmp_SampleID}.bam -c ${tmp_SampleID_trunc}-INPUT.bam -o ROSE_super_enhancer -t 2500" >> ROSE_command.sh
      chmod +x ROSE_command.sh
      # [ ! -z ${sample_H3K27AC_bam} ] && [ ! -z ${sample_INPUT_bam} ] && python ./ROSE_main.py -g HG19 -i "${tmp_sample_peak_outdir}/mergepeaks_${tmp_SampleID}.bed_goodbed.gff.better.gff" -r ${tmp_SampleID}.bam -c ${tmp_SampleID_trunc}-INPUT.bam -o ROSE_super_enhancer
      
      # run the actual ROSE command
      echo "Running the ROSE command"
      echo "Command is:"
      echo "python ./ROSE_main.py -g HG19 -i ${tmp_sample_peak_outdir}/mergepeaks_${tmp_SampleID}.bed_goodbed.gff.better.gff -r ${tmp_SampleID}.bam -c ${tmp_SampleID_trunc}-INPUT.bam -o ROSE_super_enhancer -t 2500"
      [ ! -f ROSE_super_enhancer/mergepeaks_${tmp_SampleID}_Gateway_SuperEnhancers.bed ] python ./ROSE_main.py -g HG19 -i ${tmp_sample_peak_outdir}/mergepeaks_${tmp_SampleID}.bed_goodbed.gff.better.gff -r ${tmp_SampleID}.bam -c ${tmp_SampleID_trunc}-INPUT.bam -o ROSE_super_enhancer -t 2500
      
      # try it with HOMER
      # http://homer.salk.edu/homer/ngs/peaks.html
      # findPeaks <tag directory> -i <input tag directory> -style super -o auto
      [ ! -z ${sample_H3K27AC_tagdir} ] && [ ! -z ${sample_INPUT_tagdir} ] && [ ! -f ${tmp_SampleID}_HOMER_superenhancer.txt ] && echo -e "\nRunning HOMER findPeaks:\n" && findPeaks "${sample_H3K27AC_tagdir}" -i "${sample_INPUT_tagdir}" -style super -o ${tmp_SampleID}_HOMER_superenhancer.txt -typical ${tmp_SampleID}_HOMER_typical_enhancer.txt
      # findPeaks BVI-R-H3K27AC_tagdir -i BVI-R-INPUT_tagir -style super -o HOMER_superenhancer
      # clean up the HOMER superenhancer file to be more like a real bed file
      grep -Ev '^#' ${tmp_SampleID}_HOMER_superenhancer.txt | cut -f2,3,4 > ${tmp_SampleID}_HOMER_superenhancer.bed
      grep -Ev '^#' ${tmp_SampleID}_HOMER_typical_enhancer.txt | cut -f2,3,4 > ${tmp_SampleID}_HOMER_typical_enhancer.bed
      
      
      # now overlap the new HOMER superenhancer peaks with the peaks we already have.. 
      # rm -rf HOMER_superenhancer_mergepeaks; 
      mkdir -p HOMER_superenhancer_mergepeaks
      cd HOMER_superenhancer_mergepeaks
      ln -s ../peaks/${tmp_SampleID}.bed
      ln -s ../peaks/${tmp_SampleID_H3K4ME3}.bed
      ln -s ../peaks/Gencode.bed
      ln -s ../${tmp_SampleID}_HOMER_superenhancer.txt
      [ ! -f venn.txt ] && mergePeaks ${tmp_SampleID}_HOMER_superenhancer.txt ${tmp_SampleID}.bed ${tmp_SampleID_H3K4ME3}.bed Gencode.bed -prefix mergepeaks -venn venn.txt
      [ ! -f ${tmp_SampleID}_venn.pdf ] && module switch r r/3.2.0 && Rscript /ifs/home/kellys04/projects/myBioinformatics/HOMER_mergePeaks_multiVenn/multi_peaks_Venn.R "${tmp_SampleID}" venn.txt
      [ ! -f ${tmp_SampleID}_UpSetR_plot.pdf ] && module switch r r/3.3.0 && Rscript /ifs/home/kellys04/projects/myBioinformatics/HOMER_mergePeaks_venn_UpSetR/multi_peaks_UpSet_plot.R "${tmp_SampleID}" venn.txt
      
      
      
      cd "$tmp_sample_outdir"
      
      # make a dir for all the peaks for comparison
      # echo "Making a dir to hold copies of all peaks for comparisons, dir is:"
      echo "Copying peaks to common dir:"
      echo "$common_branch_peaks"
      # mkdir -p ${branch}_${tmp_SampleID}_all_peaks_compare
      # cd ${branch}_${tmp_SampleID}_all_peaks_compare
      # echo "new pwd is:"
      # pwd
      /bin/cp ${tmp_SampleID}_HOMER_superenhancer.bed ${common_branch_peaks}/
      /bin/cp ${tmp_SampleID}_HOMER_typical_enhancer.bed ${common_branch_peaks}/
      /bin/cp ${tmp_sample_peak_outdir}/${tmp_SampleID}.bed ${common_branch_peaks}/
      /bin/cp ${tmp_sample_peak_outdir}/${tmp_SampleID_H3K4ME3}.bed ${common_branch_peaks}/
      /bin/cp -n ${tmp_sample_peak_outdir}/Gencode.bed ${common_branch_peaks}/
      
      echo -e "\npwd is:\n"
      pwd
      echo -e "\npwd contents:\n"
      ls -l *
      
      # HOMER superenhancer post processing
      # copy the HOMER super enhancers to the external dir
      # they should all be the same because the HOMER ones were called from the alignment not the MACS2 branch
      # mkdir -p ${external_super_dir}/HOMER/${branch}
      # only copy for the main branch thats all that matters for HOMER
      echo "HOMER external outdir is:"
      echo "$(readlink -m ${external_super_dir})"
      echo "${external_super_dir}/HOMER"
      mkdir -p ${external_super_dir}/HOMER
      mkdir -p ${external_super_dir}/HOMER/bed
      # make a symlink for easy access
      ln -s ${external_super_dir} ${main_outdir}/external_superenhancer_dir
      # /ifs/data/sequence/results/smithlab/ChIPSeq/2016-01-04/ChIP-Seq/superenhancer
      
      # copy just for that branch, also make a BigBed and get the URL
      # https://genome.ucsc.edu/goldenPath/help/bigBed.html
      # need to also sort the bed 
      [ "${branch}" = "peaks.by_sample.macs_broad" ] && "$bedsort_path" "${tmp_SampleID}_HOMER_superenhancer.bed" "${tmp_SampleID}_HOMER_superenhancer.bed.sorted"
      [ -f "${tmp_SampleID}_HOMER_superenhancer.bed.sorted" ] && "$bed_to_bigbed_path" "${tmp_SampleID}_HOMER_superenhancer.bed.sorted" "$hg19_chrom_sizes" "${tmp_SampleID}_HOMER_superenhancer.bigbed"
      [ -f ${tmp_SampleID}_HOMER_superenhancer.bigbed ] && /bin/cp ${tmp_SampleID}_HOMER_superenhancer.bigbed ${external_super_dir}/HOMER/${tmp_SampleID}_HOMER_superenhancer.bigbed && echo "${myurlbase}${external_super_dir##/ifs/data/sequence/results/}/HOMER/${tmp_SampleID}_HOMER_superenhancer.bigbed" > ${tmp_SampleID}_HOMER_superenhancer_url.txt
      [ "${branch}" = "peaks.by_sample.macs_broad" ] && [[ "${tmp_SampleID}" == *H3K27AC* ]] && /bin/cp "${tmp_SampleID}_HOMER_superenhancer.bed.sorted" "${external_super_dir}/HOMER/bed/${tmp_SampleID}_HOMER_superenhancer.bed.sorted"
      
      # ROSE superenhancer post processing
      # # copy the ROSE super enhancers to 
      echo -e "ROSE external dir is:"
      echo "$(readlink -m ${external_super_dir})"
      echo "${external_super_dir}/ROSE"
      mkdir -p "${external_super_dir}/ROSE"
      mkdir -p "${external_super_dir}/ROSE/bed"
      # # need to sort the bed & convert to BigBed
      [ "${branch}" = "peaks.by_sample.macs_broad" ] && [ -f ROSE_super_enhancer/mergepeaks_${tmp_SampleID}_Gateway_SuperEnhancers.bed ] && "$bedsort_path" "ROSE_super_enhancer/mergepeaks_${tmp_SampleID}_Gateway_SuperEnhancers.bed" "ROSE_super_enhancer/mergepeaks_${tmp_SampleID}_Gateway_SuperEnhancers.bed.sorted"
      [ -f ROSE_super_enhancer/mergepeaks_${tmp_SampleID}_Gateway_SuperEnhancers.bed.sorted ] && "$bed_to_bigbed_path" "ROSE_super_enhancer/mergepeaks_${tmp_SampleID}_Gateway_SuperEnhancers.bed.sorted" "$hg19_chrom_sizes" "ROSE_super_enhancer/mergepeaks_${tmp_SampleID}_Gateway_SuperEnhancers.bed.sorted.bigbed"
      [ -f ROSE_super_enhancer/mergepeaks_${tmp_SampleID}_Gateway_SuperEnhancers.bed.sorted.bigbed ] && /bin/cp ROSE_super_enhancer/mergepeaks_${tmp_SampleID}_Gateway_SuperEnhancers.bed.sorted.bigbed ${external_super_dir}/ROSE/${tmp_SampleID}_ROSE_superenhancer.bigbed && echo "${myurlbase}${external_super_dir##/ifs/data/sequence/results/}/ROSE/${tmp_SampleID}_ROSE_superenhancer.bigbed" > ${tmp_SampleID}_ROSE_superenhancer_url.txt
      [ "${branch}" = "peaks.by_sample.macs_broad" ] && [[ "${tmp_SampleID}" == *H3K27AC* ]] && /bin/cp "ROSE_super_enhancer/mergepeaks_${tmp_SampleID}_Gateway_SuperEnhancers.bed.sorted" "${external_super_dir}/ROSE/bed/${tmp_SampleID}_ROSE_superenhancer.bed"
      
      # # print the UCSC custom track for both HOMER and ROSE
      # print the ROSE custom track
      [ -f ROSE_super_enhancer/mergepeaks_${tmp_SampleID}_Gateway_SuperEnhancers.bed.sorted.bigbed ] && echo -e "${tmp_SampleID}_a\ttrack type=bigBed name=\"${tmp_SampleID}_RS:SE\" description=\"${tmp_SampleID} ROSE superenhancer\" bigDataUrl=${myurlbase}${external_super_dir##/ifs/data/sequence/results/}/ROSE/${tmp_SampleID}_ROSE_superenhancer.bigbed" > ${tmp_SampleID}_all_superenhancer_UCSC.txt
      # print the HOMER custom track
      [ -f ${tmp_SampleID}_HOMER_superenhancer.bigbed ] && echo -e "${tmp_SampleID}_b\ttrack type=bigBed name=\"${tmp_SampleID}_HM:SE\" description=\"${tmp_SampleID} HOMER superenhancer\" bigDataUrl=${myurlbase}${external_super_dir##/ifs/data/sequence/results/}/HOMER/${tmp_SampleID}_HOMER_superenhancer.bigbed" >> ${tmp_SampleID}_all_superenhancer_UCSC.txt
      # print the bigWig custom track
      [ -f ${tmp_SampleID}_HOMER_superenhancer.bigbed ] && [ -n ${sample_bigWig} ] && echo -e "${tmp_SampleID}_c\ttrack type=bigWig name=\"${tmp_SampleID}_bW\" bigDataUrl=${myurlbase}${sample_bigWig##/ifs/data/sequence/results/} ${UCSC_bigwig_view}" >>  ${tmp_SampleID}_all_superenhancer_UCSC.txt
      # add the custom tracks to the overall UCSC custom track file
      [ -f ${tmp_SampleID}_HOMER_superenhancer.bigbed ] && [ -n ${sample_bigWig} ] && cat ${tmp_SampleID}_all_superenhancer_UCSC.txt >> "${main_outdir}/UCSC_custom_track.txt.tmp"
      
      
      
E0F

    
  # ) & done
  done
  echo -e "\n\n"
  
done

# copy all the compare peaks dirs
cd $main_outdir

# clean up the UCSC sesion doc
# add the default browser postion and sort the entries, remove the first column which was used as a sort key
echo "browser position chr1:2,439,912-2,528,361" > "${main_outdir}/UCSC_custom_track.txt"
cat "${main_outdir}/UCSC_custom_track.txt.tmp" | sort | cut -f2- >> "${main_outdir}/UCSC_custom_track.txt"

# rm -rf peaks_for_browser; 
# mkdir -p peaks_for_browser
# find . -type d -name "*_all_peaks_compare*" -exec mv {} peaks_for_browser/ \;
# THE END
exit; exit; exit
# ^ in case I broke a loop or something, please please exit







# ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ 
# ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ 
# ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ 
# ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ 
# notes and stuff below:

# ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ 
# for submitting this pipeline script:
pipeline_script="/ifs/home/kellys04/projects/SmithLab_ChIPSeq_2016-03-31/project_notes/super_enhancer_calling/superenhancer_pipeline.sh"
tmp_logdir="/ifs/home/kellys04/projects/SmithLab_ChIPSeq_2016-03-31/project_notes/super_enhancer_calling/ROSE_superenhancers_2/logs"
mkdir -p "$tmp_logdir"
# qsub -o :${tmp_logdir}/ -e :${tmp_logdir}/ -pe threaded 1 -l h_vmem=10G -l mem_free=10G -l mem_token=10G "$pipeline_script"
qsub -q all.q@node039.cm.cluster -o :${tmp_logdir}/ -e :${tmp_logdir}/ -pe threaded 1 "$pipeline_script"
# ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ 
# usethe same node each time runs faster I think?
