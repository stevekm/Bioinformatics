#!/bin/bash

# USAGE: VarScan_SomaticSniper_vcf_annotate.sh
# DESCRIPTION: This script will annotate Somatic Sniper produced VCF files using ANNOVAR

# this is depricated; I've hardcoded the args for now 
# InputFile="$1"
# OutDir="$2"
# Genome="$3"


Genome="hg19"


# echo "$InputFile is InputFile"
# echo "$OutDir is OutDir"
# echo "$Genome is Genome"


### ~~~~ Settings ~~~~~~
# location of the annovar convert program
ANNOVAR_CONVERT=~/software/annovar/convert2annovar.pl

# table command
ANNOVAR_TABlE=~/software/annovar/table_annovar.pl

# location of the ANNOVAR databases; I changed the name of the dir btw
DB_DIR=~/software/annovar/db

# Convert the VCF to TSV, using libvcf <https://github.com/ekg/vcflib>
VCF_to_TSV=~/software/vcflib/bin/vcf2tsv

# location of the ANNOVAR conversion script to produce REDCap output
ANNOVAR2REDCAP=~/projects/clinical_genomic_reporting/code/annovar2redcap.R
# ANNOVAR2REDCAP=/ifs/home/$(whoami)/projects/clinical_genomic_reporting/code/annovar2redcap_old.R # use the old one while I develop the new one


ProjDir="~/projects/clinical_genomic_reporting"
Mutations_dir="$ProjDir/556_panel_reporting/NYU556/consolidated.mutations.indels"
MyScript="~/projects/clinical_genomic_reporting/556_panel_reporting/code/556_panel_ANNOVAR.sh"
OutDir="~/projects/clinical_genomic_reporting/556_panel_reporting/annotations"
Genome="hg19"

SomaticSniper_dir="~/projects/clinical_genomic_reporting/556_panel_reporting/NYU556/mutations/SS.SI"
##### ~~~~~~~~~

for i in $SomaticSniper_dir/*; do
  if [[ -d $i ]]; then # make sure its a directory

    # get sample name from dirname
    SampleName=$(basename $i)

    # get the VCFs for the sample
    SampleVCF=$i/*.vcf

    for q in $SampleVCF; do
      # echo $q
      # make the outdir for the sample
      mkdir -p $OutDir/"$SampleName"
      # convert to annovar format
      ANNOVAR_file=$OutDir/"$SampleName"/"$(basename "$q")".avinput
      if [ ! -e "$ANNOVAR_file" ]; then 
        $ANNOVAR_CONVERT -format vcf4old "$q" -includeinfo > $ANNOVAR_file
      fi
      
      # annotate
      ANNOVAR_output=${ANNOVAR_file//.vcf*/}
      if [ ! -e "${ANNOVAR_output}.${Genome}_multianno.csv" ]; then 
  	    $ANNOVAR_TABlE "$ANNOVAR_file" $DB_DIR -buildver $Genome -out "$ANNOVAR_output" -remove -protocol refGene,ljb26_all,cosmic68,clinvar_20150629, -operation g,f,f,f -nastring . -csvout
  	  fi
  
      # also convert VCF to TSV
      CONVERTED_TSV=$OutDir/"$SampleName"/"$(basename "$q")".tsv
      if [ ! -e "$CONVERTED_TSV" ]; then 
        $VCF_to_TSV "$q" > $CONVERTED_TSV
      fi
      
      
    done
  fi
done

