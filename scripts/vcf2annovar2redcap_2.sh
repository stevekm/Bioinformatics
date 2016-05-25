#!/bin/bash

##
## USAGE: vcf2annovar2redcap_2.sh /path/to/vcf_input_dir /path/to/outdir /path/to/params
## OR
##  vcf2annovar2redcap.sh /path/to/vcf_file.vcf /path/to/outdir /path/to/params
## annotate all vcf's within the directory, output a table in the outdir
## developer's note: never ever use . in dirnames you you'll probably break this script

# This script will:
# convert VCF files in the directory into a format that can be used by ANNOVAR for annotation
# then, annotate the file with ANNOVAR
# then, convert it to a TSV because thats handy
# then, parse the ANNOVAR output in R and write out a CSV in the format we need for REDCap
# output all files in the same directory containing VCFs

# turn on extglob to find files better
# shopt -s extglob
# set -B
# set -x

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
  echo "PARAMS NOT FOUND, CHECK FILE PATH"
fi

# make the outdir if it doesn't exist
mkdir -p $OUT_DIR

# adapt the script to use either a file or dir
# # check if its a file
if [ -f $1 ]; then
  echo "ITS A FILE"
  VCF_DIR=$(dirname $1) # set the parent dir
  INPUT_ITEM=$1
# if its a dir instead..
elif [ -d $1 ]; then
  echo "Its A DIR"
  VCF_DIR=$1
  INPUT_ITEM="*.vcf" # set the search string for the .vcf's
  
  # check to make sure dir contains the VCF's
	# # detect if only .vcf.gz files are present and unzip them
	# # # look for absence of .vcf, presence of .vcf.gz
	if [[ -z $(find $VCF_DIR -type f -name "*.vcf") && -n $(find $VCF_DIR -type f -name "*.vcf.gz") ]]; then
	  echo "NO VCF FOUND, GZ FOUND"
	  echo "UNZIPPING GZ"
	  for i in $(find $VCF_DIR -type f -name "*.vcf.gz"); do
	    # strip .gz from name
	    zz=${i//.gz/}
	    gunzip -vc $i > $zz
	  done
	# look for absence of VCF only
	elif [[ -z $(find $VCF_DIR -type f -name "*.vcf") ]]; then
	  echo "NO VCF's FOUND"
	# look for presence of .vcf
	elif [[ -n $(find $VCF_DIR -type f -name "*.vcf") ]]; then
	  echo "VCF FOUND"
	fi
fi

echo -e "Dir is:\n$VCF_DIR"
echo "Item is $INPUT_ITEM"
echo -e "Command is:" #\nfind $VCF_DIR -path ${INPUT_ITEM} -type f"
echo -e "find $VCF_DIR -path "$INPUT_ITEM" -type f"

# Annotate the VCF's with ANNOVAR
for i in $(find $VCF_DIR -path "$INPUT_ITEM" -type f); do 
  echo -e "Processing file $i\n"

	# strip the file extension from the name 
	zz=${i//.vcf/}

	# convert the vcf to annovar input format
	if [[ -f ${i}.avinput ]]; then
	  echo "ANNOVAR INPUT ALREADY EXISTS, SKIPPING CONVERSION STEP"
	else
	  echo "Converting $(basename $i) to ANNOVAR format"
	  echo "${ANNOVAR_CONVERT} -format vcf4old ${i} -includeinfo > ${i}.avinput"
	  ${ANNOVAR_CONVERT} -format vcf4old ${i} -includeinfo > ${i}.avinput
	fi

	# annotate with ANNOVAR table
	if [[ -f $zz.${BUILD_VERSION}_multianno.csv ]]; then
	  echo "ANNOVAR TABLE ALREADY EXISTS, SKIPPING STEP"
	else
	  echo "Annotating $(basename ${i}.avinput) with ANNOVAR"
	  echo "${ANNOVAR_TABlE} ${i}.avinput ${DB_DIR} -buildver $BUILD_VERSION -out ${zz} -remove -protocol refGene,ljb26_all,cosmic68 -operation g,f,f -nastring . -csvout"
	  ${ANNOVAR_TABlE} ${i}.avinput ${DB_DIR} -buildver $BUILD_VERSION -out ${zz} -remove -protocol refGene,ljb26_all,cosmic68 -operation g,f,f -nastring . -csvout
	fi

	# convert the VCF to TSV for later
	if [[ -f ${i}.tsv ]]; then
	  echo "TSV ALREADY EXISTS, SKIPPING STEP"
	else
	  echo "Converting to TSV"
	  echo "${VCF_to_TSV} ${i} > ${i}.tsv"
	  ${VCF_to_TSV} ${i} > ${i}.tsv
	fi
  
  # pass to R for processing for REDCap.. 
  echo "Making table to upload to REDCap"
  echo "${ANNOVAR2REDCAP} $(dirname $i)/$(basename $zz).${BUILD_VERSION}_multianno.csv"
	${ANNOVAR2REDCAP} $(dirname $i)/$(basename $zz).${BUILD_VERSION}_multianno.csv
done
