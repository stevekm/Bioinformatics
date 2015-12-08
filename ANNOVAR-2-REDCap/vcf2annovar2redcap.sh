#!/bin/bash

##
## USAGE: vcf2annovar2redcap.sh /path/to/vcf_input_dir /path/to/outdir /path/to/params
## OR
##  vcf2annovar2redcap.sh /path/to/vcf_file.vcf /path/to/outdir /path/to/params
## annotate all vcf's within the directory, output a table in the outdir
## Input can be a single .vcf or a dir with many of these, or a zip containing vcf or vcf.gz
## .vcf.gz will automatically be extracted before processing
## EXAMPLE:
## ./code/vcf2annovar2redcap.sh ./sample_official_VCF ./redcap_project_files ./code/vcf2annovar2redcap_params.sh

# This script will:
# convert VCF files in the directory into a format that can be used by ANNOVAR for annotation
# then, annotate the file with ANNOVAR
# then, convert it to a TSV 
# then, parse the ANNOVAR output in R and write out a CSV in the format we need for REDCap
# output all files in the same directory containing VCFs

# process command-line inputs
# # if not enough args, output USAGE and exit
if (($# != 3)); then
  grep '^##' $0
  exit
fi

# # Input Script Args
VCF_DIR=$1 # depricated, see below.. 
OUT_DIR=$2
PARAMS=$3

# get the params from file (required!)
if [ -f $PARAMS ]; then 
  source $PARAMS
else
  echo -e "\nParams not found, try checking the file path"
fi

# make the outdir if it doesn't exist
mkdir -p $OUT_DIR

# set the name for the final output table
OUTPUT_GENE_TABLE=${OUT_DIR}/gene_table_$(date -u +%Y%m%dt%H%M%S).csv

# #
# # INPUT FILE CHECKING
# #
# # check if its a file
if [ -f $1 ]; then
  echo -e "\nFile found"
  VCF_DIR=$(dirname $1) # set the parent dir
  INPUT_ITEM=$1
  
  # check the file extension; make sure its a VCF, if its a .vcf.gz then unzip it
  # # #
  # # #
  # # # look for absence of .vcf, presence of .vcf.gz 
  if [[ $INPUT_ITEM != *.vcf && $INPUT_ITEM == *.vcf.gz ]]; then
    echo -e "\nFile is not a VCF\nFile is a VCF.GZ; Uzipping"
    # strip the file extension
    zz=${INPUT_ITEM%%.gz} 
    gunzip -vc $INPUT_ITEM > $zz
    # set the $INPUT_ITEM to zz
    INPUT_ITEM=$zz
    # find the files to annotate
    FILES=$(find $VCF_DIR -path "$INPUT_ITEM" -type f)
    
  # # #
  # # #
  # # # look for absence of .vcf, presence of .vcf.zip
  elif [[ $INPUT_ITEM != *.vcf && $INPUT_ITEM == *.vcf.zip ]]; then
    echo -e "\nFile is not a VCF\nFile is a VCF.ZIP; Unzipping"
    unzip -o $INPUT_ITEM -d $(dirname $INPUT_ITEM)
    # find the files to annotate
    # #  look for presence of vcf's
    if [[ -n $(find $VCF_DIR -iname "*.vcf" -type f) ]]; then
      echo -e "\nUnzipped VCF files found, Proceeding with annotation"
      # set the list of vcf's
      FILES=$(find $VCF_DIR -iname "*.vcf" -type f)
    # # look for absence of vcf's, presence of vcf.gz's
    elif [[ -z $(find $VCF_DIR -iname "*.vcf" -type f) && -n $(find $VCF_DIR -iname "*.vcf.gz" -type f) ]]; then
      echo -e "\nVCF files not found\nVCF.GZ files found; Unzipping"
      gunzip -vf $VCF_DIR/*.vcf.gz
      # set the list of vcf's
      FILES=$(find $VCF_DIR -iname "*.vcf" -type f)
    # # if no vcf or vcf.gz found  
    else
      echo -e "\nNeither VCF nor VCF.GZ files found\nExiting"
      exit
    fi
    
  # # #
  # # #
  # # # look for absence of .vcf only
  elif [[ $INPUT_ITEM != *.vcf ]]; then
    echo -e "\nFile is not a VCF\nExiting"
    exit
  
  # # #
  # # #
  # # # look for presence of .vcf
  elif [[ $INPUT_ITEM == *.vcf ]]; then
    echo -e "\nVCF file found\nProceeding with annotation"
    
  else
    echo -e "\nFile not recognized\nExiting"
    exit
  fi

# #  
# # if its a dir instead..
elif [ -d $1 ]; then
  echo -e "\nDirectory found"
  VCF_DIR=$1
  INPUT_ITEM="*.vcf" # set the search string for the .vcf's # this is depricated in this case

  # # #
  # # #
  # # # look for .gz, unzip them
  for i in $(find $VCF_DIR -type f -name "*.vcf.gz"); do
    # strip .gz from name
	  zz=${i%%.gz} 
	  gunzip -vc $i > $zz
  done

#   if [[ -z $(find $VCF_DIR -type f -name "*.vcf") ]]; then
#     echo -e "\nVCF files not found\nExiting"
#     exit
#   fi
# 

  # get the VCF files
  FILES=$(find $VCF_DIR -type f -name "*.vcf")
  
  # # #
  # # #
  # # # make sure there really are VCF files present
  if [[ -z $FILES ]]; then
    echo -e "\nVCF files not found\nExiting"
    exit
  fi
  
fi

echo -e "\nInput directory is:\n$VCF_DIR"
# echo -e "\nITEM IS:\n $INPUT_ITEM"
# echo -e "\nCOMMAND IS:"
# echo -e "\nfind $VCF_DIR -path "$INPUT_ITEM" -type f"
echo -e "\nFiles found:"
# find $VCF_DIR -path "$INPUT_ITEM" -type f

for i in $FILES; do
  echo $i
done

#
# # 
# # # 
# Annotate the VCF's with ANNOVAR
# for i in $(find $VCF_DIR -path "$INPUT_ITEM" -type f); do 
for i in $FILES; do 
  
  echo -e "\n---------------\n"
  echo -e "\nProcessing file $i\n"

	# strip the file extension from the name 
	zz=${i//.vcf/}

	# convert the vcf to annovar input format
	if [[ -f ${i}.avinput ]]; then
	  echo -e "\nANNOVAR input already exists, skipping conversion step"
	else
	  echo -e "\nConverting $(basename $i) to ANNOVAR format\nCommand:"
	  echo -e "\n${ANNOVAR_CONVERT} -format vcf4old ${i} -includeinfo > ${i}.avinput"
	  ${ANNOVAR_CONVERT} -format vcf4old ${i} -includeinfo > ${i}.avinput
	fi

	# annotate with ANNOVAR table
	if [[ -f $zz.${BUILD_VERSION}_multianno.csv ]]; then
	  echo -e "\nANNOVAR table already exists, skipping annotation step"
	else
	  echo -e "\nAnnotating $(basename ${i}.avinput) with ANNOVAR\nCommand:"
	  echo -e "\n${ANNOVAR_TABlE} ${i}.avinput ${DB_DIR} -buildver $BUILD_VERSION -out ${zz} -remove -protocol refGene,ljb26_all,cosmic68 -operation g,f,f -nastring . -csvout"
	  ${ANNOVAR_TABlE} ${i}.avinput ${DB_DIR} -buildver $BUILD_VERSION -out ${zz} -remove -protocol refGene,ljb26_all,cosmic68 -operation g,f,f -nastring . -csvout
	fi

	# convert the VCF to TSV for later
	if [[ -f ${i}.tsv ]]; then
	  echo -e "\nTSV file already exists for this sample, skipping conversion step"
    # # get the TSV file
    CONVERTED_TSV=${i}.tsv
    echo -e "\nTSV file is $CONVERTED_TSV"
	else
	  echo -e "\nConverting VCF to TSV"
	  CONVERTED_TSV=${i}.tsv
	  echo -e "\nCommand is:\n"
	  echo -e "\n${VCF_to_TSV} ${i} > $CONVERTED_TSV"
	  ${VCF_to_TSV} ${i} > $CONVERTED_TSV
	fi
  
  # pass to R for processing for REDCap.. 
  # # get the TSV file
  CONVERTED_TSV=${i}.tsv

  # # get the file name
  ANNOTATED_CSV=$(dirname $i)/$(basename $zz).${BUILD_VERSION}_multianno.csv
  
  # # get the ANNOVAR version
  ANNOVAR_VERSION=$( $ANNOVAR_TABlE | grep Version | tr -d "," | tr -d '$' | sed 's/^ *Version: Date: *//' | xargs )

  # create a column that contains the name of the genetable.csv file
  OUTPUT_GENE_TABLE_BASENAME=$(basename $OUTPUT_GENE_TABLE)

  
  # # pass the items to the R script
  # # # check if $1 was a zip; if so, pass it as well as the run name
  if [[ $1 == *.zip ]]; then
    RUN_NAME=$(basename $1)
    RUN_NAME=${RUN_NAME%%.*}
  else
    RUN_NAME="."
  fi

  echo -e "\nMaking gene table for upload to REDCap"
  echo -e "\nInput file is:\n${1}"
  echo -e "\nRun name is:\n${RUN_NAME}"
  echo -e "\nR command is:\n${ANNOVAR2REDCAP} $ANNOTATED_CSV $i "$ANNOVAR_VERSION" $CONVERTED_TSV $OUTPUT_GENE_TABLE_BASENAME "$RUN_NAME""
	${ANNOVAR2REDCAP} $ANNOTATED_CSV $i "$ANNOVAR_VERSION" $CONVERTED_TSV $OUTPUT_GENE_TABLE_BASENAME "$RUN_NAME"
	
done


#
# # 
# # # 
# concatenate all the gene_table.csv files produced
echo -e "\nWriting gene table to file:\n${OUTPUT_GENE_TABLE}\n"

q=0                                       # Reset a counter
for i in $VCF_DIR/*_gene_table.csv; do 
  # echo $i
   if [[ $q -eq 0 ]] ; then 
      head -1  $i >   $OUTPUT_GENE_TABLE # Copy header if it is the first file
   fi
   tail -n +2  $i >>  $OUTPUT_GENE_TABLE # Append from the 2nd line each file
   q=$(( $q + 1 ))                        # Increase the counter
done


#
# # 
# # # 
# upload the file # code works, but server has issues
echo -e "\nUploading the gene table file\n"
OUTPUT_GENE_TABLE_fullpath=$(readlink -m $OUTPUT_GENE_TABLE)
echo -e "\n${REDCAP_IMPORT_R} ${OUTPUT_GENE_TABLE_fullpath} SERVICE TOKEN\n"
${REDCAP_IMPORT_R} ${OUTPUT_GENE_TABLE_fullpath} "${SERVICE}" ${TOKEN}


