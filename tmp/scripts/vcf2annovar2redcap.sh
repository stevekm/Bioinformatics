#!/bin/bash
# USAGE: vcf2annovar2redcap.sh /path/to/vcf_dir
# convert VCF files in the directory into a format that can be used by ANNOVAR for annotation
# then, annotate the file with ANNOVAR
# then, convert it to a TSV because thats handy
# then, parse the ANNOVAR output in R and write out a CSV in the format we need for REDCap
# output all files in the same directory containing VCFs

# PATH to programs used
# location of the annovar convert program
ANNOVAR_CONVERT=/ifs/home/username/software/annovar/convert2annovar.pl

# table command
ANNOVAR_TABlE=/ifs/home/username/software/annovar/table_annovar.pl

# location of the databases; I changed the name of the dir btw
DB_DIR=/ifs/home/username/software/annovar/db

# Convert the VCF to TSV, using libvcf <https://github.com/ekg/vcflib>
VCF_to_TSV=~/software/vcflib/bin/vcf2tsv

ANNOVAR2REDCAP=/ifs/home/username/projects/clinical_genomic_reporting/code/annovar2redcap.R
# location of the test VCF to start with

# get the VCF dir
VCF_DIR=$1

cd $VCF_DIR
for i in $(find $VCF_DIR -type f -name "*.vcf"); do 
# $i # include full file path

# echo $i
# strip the file extension for later
zz=$(basename $i)
zz=${zz//.vcf/}
# echo $zz.
# head $i

# convert the vcf to annovar input format
if [[ -f ${i}.avinput ]]; then
  echo "ANNOVAR INPUT ALREADY EXISTS, SKIPPING STEP"
else
  ${ANNOVAR_CONVERT} -format vcf4old ${i} -includeinfo > ${i}.avinput
fi


# annotate with ANNOVAR table
if [[ -f $(dirname $i)/$zz.hg19_multianno.csv ]]; then
  echo "ANNOVAR TABLE ALREADY EXISTS, SKIPPING STEP"
else
  ${ANNOVAR_TABlE} ${i}.avinput ${DB_DIR} -buildver hg19 -out ${zz} -remove -protocol refGene,ljb26_all,cosmic68 -operation g,f,f -nastring . -csvout
fi


# convert the VCF to TSV for later
if [[ -f ${i}.tsv ]]; then
echo "TSV ALREADY EXISTS, SKIPPING STEP"
else
  ${VCF_to_TSV} ${i} > ${i}.tsv
  # echo ${ANNOVAR2REDCAP} $(dirname $i)/$zz.hg19_multianno.csv
fi

${ANNOVAR2REDCAP} $(dirname $i)/$zz.hg19_multianno.csv

done