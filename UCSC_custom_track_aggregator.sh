#!/bin/bash
set -x

## USAGE: UCSC_custom_track_aggregator.sh
## this script will aggregate all bed files for use in UCSC, convert them to bigBed format, and
## create a table with UCSC custom browser track entries for each file per sample
## warning: find steps depend on the presence of only one matching file; adjust accordingly !!

# dir to hold all the intermediary files
ProjDir="/ifs/home/kellys04/projects/SmithLab_ChIP-Seq_2016-03-31/project_notes/UCSC_custom_tracks"
mkdir -p "$ProjDir"
cd "$ProjDir"

tmp_track_dir="${ProjDir}/tmp_tracks"
# remove any files present
find ${tmp_track_dir} -type f -exec rm -f {} \;
mkdir -p "$tmp_track_dir"

# our Gencode TSS regions bed filea; TSS regions +/-10kbp
gen_bed="/ifs/home/kellys04/projects/SmithLab_ChIP-Seq_2016-03-31/project_data2/gencode.v19.annotation_TSSd500_10kbp.bed"

# dir with the peak calling results; one subdir per sample, sample ID  = subdir name
pipeline_peaks_dir="/ifs/home/kellys04/projects/SmithLab_ChIP-Seq_2016-03-31/pipeline/peaks/results/peaks.by_sample.macs_broad/align.by_sample.bowtie2"

# get a list of the sample ID's from the pipeline peaks results dir
sampleIDs=$(find "$pipeline_peaks_dir" -maxdepth 1 -mindepth 1 -type d -exec basename {} \;  | sort -u) # AGK-D-H3K27AC ...
patientIDs=$(find "$pipeline_peaks_dir" -maxdepth 1 -mindepth 1 -type d -exec basename {} \;  | sed 's/-.*//g' | sort -u) # AGK BVI CBK DKJ EVJ FLV GHW IDY PRW SPN SRR ZGR ZNK
markIDs=$(find "$pipeline_peaks_dir" -maxdepth 1 -mindepth 1 -type d -exec basename {} \;  | sed 's/^\(.*-\)\([HI].*\)$/\2/g' | sort -u) # H3K27AC H3K27ME3 H3K4ME3 H3K9AC H3K9ME3 INPUT
declare -a sampleIDs_array; declare -a patientIDs_array; declare -a markIDs_array
sampleIDs_array=($sampleIDs)
patientIDs_array=($patientIDs)
markIDs_array=($markIDs)

# file to hold the custom tracks for UCSC
UCSC_custom_track_file="${ProjDir}/UCSC_custom_track.txt"
# clear the UCSC custom track file
echo -n "" > ${UCSC_custom_track_file}.tmp
# starting position for the browser
browser_pos='browser position chr1:2,439,912-2,528,361'

# file path to directory containing the external files for use in the browser
# # these already exist and contain our files !!
external_results_dir="/ifs/data/sequence/results/smithlab/ChIP-Seq"; mkdir -p "$external_results_dir"
external_peaks_dir="${external_results_dir}/peaks"; mkdir -p "$external_peaks_dir"
external_super_dir="${external_results_dir}/superenhancer"; mkdir -p "$external_super_dir"
external_bigwig_dir="${external_results_dir}/alignment/bigwigs"; mkdir -p "$external_bigwig_dir"


# symlink to the external dir from the proj dir
ln -fs "$external_results_dir" "${ProjDir}/external_results_dir"
# some items & programs for file conversions
hg19_chrom_sizes="/ifs/home/kellys04/software/hg19.chrom.sizes"
bed_to_bigbed="/ifs/home/kellys04/software/UCSC/bedToBigBed"
bedsort="/ifs/home/kellys04/software/UCSC/bedSort"


# ~~~~~~~ UCSC CUSTOM TRACK FUNCTIONS TO USE ~~~~~~~~~~~ # 

get_external_URL()
{
  ## USAGE: get_external_URL file
  ## file must be in the external dir given below!
  func_resultspath_base="/ifs/data/sequence/results/"
  func_urlbase="https://genome.med.nyu.edu/results"
  func_path=$(readlink -m "$1")
  echo "${func_urlbase}/${func_path##$func_resultspath_base}"
}

ucsc_cust_track_bigBed()
{
  ## USAGE: ucsc_cust_track_bigBed 'track name' 'track description' 'track URL'
  func_tracktype='bigBed'
  func_trackname="$1"
  func_trackdescr="$2"
  func_trackURL="$3"
  bigbed_track_template="track type=${func_tracktype} name=\"${func_trackname}\" description=\"${func_trackdescr}\" bigDataUrl=${func_trackURL}"
  echo "$bigbed_track_template"
}

ucsc_cust_track_bigWig()
{
  ## USAGE: ucsc_cust_track_bigWig 'track name' 'track URL'
  func_tracktype="bigWig"
  func_trackname="$1"
  func_trackURL="$2"
  func_bigwig_view='visibility=full autoScale=off alwaysZero=on maxHeightPixels=50 graphType=bar viewLimits=0:0.3'
  bigwig_track_template="track type=${func_tracktype} name=\"${func_trackname}\" bigDataUrl=${func_trackURL} ${func_bigwig_view}"
  echo "$bigwig_track_template"
}


ROSE_track_descr()
{
  ## USAGE: ROSE_track_descr sampleID
  func_sampleID="$1"
  echo "${func_sampleID} ROSE Superenhancers"
}

HOMER_track_descr()
{
  ## USAGE HOMER_track_descr sampleID
  func_sampleID="$1"
  echo "${func_sampleID} HOMER Superenhancers"
}

transpose_awk()
{
  ## USAGE: transpose_awk file
  ## OR: cat file | transpose_awk -
  awk '
BEGIN { FS=OFS="\t" }
{
    for (rowNr=1;rowNr<=NF;rowNr++) {
        cell[rowNr,NR] = $rowNr
    }
    maxRows = (NF > maxRows ? NF : maxRows)
    maxCols = NR
}
END {
    for (rowNr=1;rowNr<=maxRows;rowNr++) {
        for (colNr=1;colNr<=maxCols;colNr++) {
            printf "%s%s", cell[rowNr,colNr], (colNr < maxCols ? OFS : ORS)
        }
    }
}
' "$1"

# output:
# X       row1    row2    row3    row4
# column1 0       3       6       9
# column2 1       4       7       10
# column3 2       5       8       11
# http://stackoverflow.com/questions/1729824/transpose-a-file-in-bash
}


# ~~~~~~~~~~~~~~~~~~ # 



# ~~~~~~~ START DOING THINGS HERE ~~~~~~~~~~~ # 
# 
# ~~~~~~~ PROCESS GENCODE TSS FILES ~~~~~~~~~~~ # 
# final output file to be made:
mkdir -p "${ProjDir}/Gencode"
gen_Bigbed="${ProjDir}/Gencode/Gencode_TSS.bigBed"
# copy over the gencode bed
[ -f $gen_bed ] && [ ! -f ${ProjDir}/Gencode/Gencode_TSS.bed ] && cat "$gen_bed" | cut -f1-3 > "${ProjDir}/Gencode/Gencode_TSS.bed"
# sort the gencode bed
[ ! -f ${ProjDir}/Gencode/Gencode_TSS.bed.sorted ] && "$bedsort" "${ProjDir}/Gencode/Gencode_TSS.bed" "${ProjDir}/Gencode/Gencode_TSS.bed.sorted"
# convert to bigBed
[ ! -f $gen_Bigbed ] && "$bed_to_bigbed" "${ProjDir}/Gencode/Gencode_TSS.bed.sorted" "$hg19_chrom_sizes" "$gen_Bigbed"
# copy to the external dir
external_gencode_TSS_bigbed="${external_peaks_dir}/Gencode_TSS.bigBed"
[ -f $gen_Bigbed ] && /bin/cp "$gen_Bigbed" "$external_gencode_TSS_bigbed" # && ucsc_cust_track_bigBed "$(basename $external_gencode_TSS_bigbed)" 'Gencode TSS regions' "$(get_external_URL "$external_gencode_TSS_bigbed")"
# ~~~~~~~~~~~~~~~~~~ # 



# ~~~~~~~ PROCESS THE H3K4ME3 and H3K27AC PEAK FILES ~~~~~~~~~~~ # 
# need to convert BED to bigBed and copy to external dir
# local dir to hold the intermediary conversion files:
mkdir -p "${ProjDir}/Peaks"
# make an external dir for the beds
mkdir -p "${external_peaks_dir}/bed"

# iterate over all of the samples; look for H3K4ME3 and H3K27AC peaks
for ((i=0; i < ${#sampleIDs_array[@]}; ++i)); do
  # check if the sample ID matches one of them
  if [[ "${sampleIDs_array[$i]}" == *H3K4ME3* ]]; then
    # get the sample ID
    tmp_sampleID="${sampleIDs_array[$i]}"
    # find the peaks file
    tmp_sample_H3K4ME3_peaks=$(find $pipeline_peaks_dir -name "peaks.bed" -path "*/${tmp_sampleID}/*")
    # if peaks file found; sort bed, keep only first 3 columns, then convert to bigBed, then copy to external dir
    [ ! -z $tmp_sample_H3K4ME3_peaks ] && "$bedsort" "$tmp_sample_H3K4ME3_peaks" "${ProjDir}/Peaks/${tmp_sampleID}_peaks.bed.sorted" && cut -f1-3 "${ProjDir}/Peaks/${tmp_sampleID}_peaks.bed.sorted" > "${ProjDir}/Peaks/${tmp_sampleID}_peaks.bed.sorted.f13" && "$bed_to_bigbed" "${ProjDir}/Peaks/${tmp_sampleID}_peaks.bed.sorted.f13" "$hg19_chrom_sizes" "${ProjDir}/Peaks/${tmp_sampleID}_peaks.bigBed"
    [ -f ${ProjDir}/Peaks/${tmp_sampleID}_peaks.bigBed ] && [ ! -f ${external_peaks_dir}/${tmp_sampleID}_peaks.bigBed ] && /bin/cp "${ProjDir}/Peaks/${tmp_sampleID}_peaks.bigBed" "${external_peaks_dir}/${tmp_sampleID}_peaks.bigBed"
    
    # also copy the bed file to the external dir 
    [ ! -f ${external_peaks_dir}/bed/${tmp_sampleID}_peaks.bed.sorted ] && /bin/cp "${ProjDir}/Peaks/${tmp_sampleID}_peaks.bed.sorted" "${external_peaks_dir}/bed/${tmp_sampleID}_peaks.bed.sorted"
    
  elif [[ "${sampleIDs_array[$i]}" == *H3K27AC* ]]; then
    # echo "Found H3K27AC; ${sampleIDs_array[$i]}"
    # get the sample ID
    tmp_sampleID="${sampleIDs_array[$i]}"
    # find the peaks file
    tmp_sample_H3K27AC_peaks=$(find $pipeline_peaks_dir -name "peaks.bed" -path "*/${tmp_sampleID}/*")
    # if peaks file found; sort, then convert to bigBed, then copy to external dir
    [ ! -z $tmp_sample_H3K27AC_peaks ] && "$bedsort" "$tmp_sample_H3K27AC_peaks" "${ProjDir}/Peaks/${tmp_sampleID}_peaks.bed.sorted" && cut -f1-3 "${ProjDir}/Peaks/${tmp_sampleID}_peaks.bed.sorted" > "${ProjDir}/Peaks/${tmp_sampleID}_peaks.bed.sorted.f13" && "$bed_to_bigbed" "${ProjDir}/Peaks/${tmp_sampleID}_peaks.bed.sorted.f13" "$hg19_chrom_sizes" "${ProjDir}/Peaks/${tmp_sampleID}_peaks.bigBed"
    [ -f ${ProjDir}/Peaks/${tmp_sampleID}_peaks.bigBed ] && [ ! -f ${external_peaks_dir}/${tmp_sampleID}_peaks.bigBed ] && /bin/cp "${ProjDir}/Peaks/${tmp_sampleID}_peaks.bigBed" "${external_peaks_dir}/${tmp_sampleID}_peaks.bigBed"
    [ ! -f ${external_peaks_dir}/bed/${tmp_sampleID}_peaks.bed.sorted ] && /bin/cp "${ProjDir}/Peaks/${tmp_sampleID}_peaks.bed.sorted" "${external_peaks_dir}/bed/${tmp_sampleID}_peaks.bed.sorted"
  fi
  
done
# ~~~~~~~~~~~~~~~~~~ # 

# ~~~~~~~ GET THE FILES FOR EACH SAMPLE ~~~~~~~~~~~ # 
# iterate over all of the samples; but only use the H3K27AC ones for tracks!

for ((i=0; i < ${#sampleIDs_array[@]}; ++i)); do
  if [[ "${sampleIDs_array[$i]}" == *H3K27AC* ]]; then
  tmp_sampleID="${sampleIDs_array[$i]}"
  # echo "${tmp_sampleID}"
  # get the HOMER and ROSE super enhancer files, and the bigWig file, and the regular peaks files
  tmp_HOMER_bigbed=$(find ${external_super_dir}/HOMER -maxdepth 1 -type f -name "*${tmp_sampleID}*")
  tmp_ROSE_bigbed=$(find ${external_super_dir}/ROSE -maxdepth 1 -type f -name "*${tmp_sampleID}*")
  tmp_bigWig=$(find ${external_bigwig_dir} -maxdepth 1 -type f -name "*${tmp_sampleID}*")
  # external_peaks_dir H3K27AC H3K4ME3 sed 's/-.*//g'
  tmp_H3K27AC_bigBed=$(find ${external_peaks_dir} -maxdepth 1 -type f -name "*H3K27AC*" -name "*$(echo $tmp_sampleID | sed 's/^\(.*-[DR]\).*$/\1/g')*")
  tmp_H3K4ME3_bigBed=$(find ${external_peaks_dir} -maxdepth 1 -type f -name "*H3K4ME3*" -name "*$(echo $tmp_sampleID | sed 's/^\(.*-[DR]\).*$/\1/g')*")
  
  # if files were found, make the custom track entries for each
  [ ! -z $tmp_HOMER_bigbed ] && tmp_HOMER_track=$(ucsc_cust_track_bigBed "$(basename $tmp_HOMER_bigbed)" "$(HOMER_track_descr $tmp_sampleID)" "$(get_external_URL $tmp_HOMER_bigbed)")
  [ ! -z $tmp_ROSE_bigbed ] && tmp_ROSE_track=$(ucsc_cust_track_bigBed "$(basename $tmp_ROSE_bigbed)" "$(ROSE_track_descr $tmp_sampleID)" "$(get_external_URL $tmp_ROSE_bigbed)")
  [ ! -z $tmp_bigWig ] && tmp_bigWig_track=$(ucsc_cust_track_bigWig "$(basename $tmp_bigWig)" "$(get_external_URL $tmp_bigWig)")
  [ ! -z $tmp_H3K27AC_bigBed ] && tmp_H3K27AC_bigBed_track=$(ucsc_cust_track_bigBed "$(basename $tmp_H3K27AC_bigBed)" "$(echo $tmp_sampleID | sed 's/^\(.*-[DR]\).*$/\1/g') H3K27AC peaks" "$(get_external_URL $tmp_H3K27AC_bigBed)")
  [ ! -z $tmp_H3K4ME3_bigBed ] && tmp_H3K4ME3_bigBed_track=$(ucsc_cust_track_bigBed "$(basename $tmp_H3K4ME3_bigBed)" "$(echo $tmp_sampleID | sed 's/^\(.*-[DR]\).*$/\1/g') H3K4ME3 peaks" "$(get_external_URL $tmp_H3K4ME3_bigBed)")
  
  # make the custom Gencode track entry as well
  [ ! -z $tmp_HOMER_bigbed ] && [ ! -z $tmp_ROSE_bigbed ] && [ -f $gen_Bigbed ] && tmp_Gencode_track=$(ucsc_cust_track_bigBed "$(basename $external_gencode_TSS_bigbed)" 'Gencode TSS regions' "$(get_external_URL "$external_gencode_TSS_bigbed")")
  
  # put all the tracks together
  [ ! -z $tmp_HOMER_bigbed ] && [ ! -z $tmp_ROSE_bigbed ] && [ -f $gen_Bigbed ] && [ ! -z $tmp_bigWig ] && [ ! -z $tmp_H3K27AC_bigBed ] && [ ! -z $tmp_H3K4ME3_bigBed ] && tmp_all_tracks=$(echo -e "${tmp_sampleID}\t${tmp_HOMER_track}\t${tmp_ROSE_track}\t${tmp_H3K27AC_bigBed_track}\t${tmp_H3K4ME3_bigBed_track}\t${tmp_Gencode_track}\t${tmp_bigWig_track}")
  # this way includes column/row labels btw..
  # [ ! -z $tmp_HOMER_bigbed ] && [ ! -z $tmp_ROSE_bigbed ] && [ -f $gen_Bigbed ] && tmp_all_tracks=$(echo -e ".\tHOMER track\tROSE track\tGencode track\tbigWig alignment\n${tmp_sampleID}\t${tmp_HOMER_track}\t${tmp_ROSE_track}\t${tmp_Gencode_track}")
  
  # write out the custom track for the sample
  [ ! -z "$tmp_all_tracks" ]  && echo "$tmp_all_tracks" > "${tmp_track_dir}/${tmp_sampleID}_tmp_track.txt" # && echo "$tmp_all_tracks"
  # [ ! -z "$tmp_all_tracks" ] && echo "$tmp_all_tracks" | transpose_awk -
  
  # echo -e "\n\n"
  fi
done

# cat all the tracks together, with the labels, and transpose
echo -en ".\tHOMER Superenhancer track\tROSE Superenhancer track\tH3K27AC peaks\tH3K4ME3 peaks\tGencode TSS region track\tbigWig alignment track\n" | cat - ${tmp_track_dir}/*_tmp_track.txt | transpose_awk - > ${ProjDir}/UCSC_custom_track.tsv

# ~~~~~~~~~~~~~~~~~~ # 


exit
# ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ 
# resources:
https://genome.ucsc.edu/
https://genome.ucsc.edu/goldenpath/help/bigBed.html
https://genome.ucsc.edu/goldenpath/help/hgTracksHelp.html
https://genome.ucsc.edu/goldenpath/help/customTrack.html

https://genome.ucsc.edu/cgi-bin/hgCustom
