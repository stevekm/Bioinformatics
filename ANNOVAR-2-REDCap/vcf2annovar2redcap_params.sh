#!/bin/bash

# params file for the vcf2annovar2redcap script
# these are relatively static configurations that need to be easy to modify without messing with the main scripts

# location of the annovar convert program
ANNOVAR_CONVERT=/ifs/home/$(whoami)/software/annovar/convert2annovar.pl

# table command
ANNOVAR_TABlE=/ifs/home/$(whoami)/software/annovar/table_annovar.pl

# location of the ANNOVAR databases; I changed the name of the dir btw
DB_DIR=/ifs/home/$(whoami)/software/annovar/db

# Convert the VCF to TSV, using libvcf <https://github.com/ekg/vcflib>
VCF_to_TSV=~/software/vcflib/bin/vcf2tsv

# location of the ANNOVAR conversion script to produce REDCap output
ANNOVAR2REDCAP=/ifs/home/$(whoami)/projects/code/annovar2redcap.R
# ANNOVAR2REDCAP=/ifs/home/$(whoami)/projects/code/annovar2redcap_old.R # use the old one while I develop the new one


# ANNOVAR build version, e.g. hg19, hg38, dm3, mm9, ... 
BUILD_VERSION=hg19

REDCap_Data_Dictionary=/ifs/home/$(whoami)/projects/data_dictionary/REDCapProject_DataDictionary_2015-12-23.csv

# Set secret token specific to your REDCap project
TOKEN=$(cat /ifs/home/$(whoami)/projects/redcap_clinical_reporting_token.txt)


# Set the url to the api (ex. https://YOUR_REDCAP_INSTALLATION/api/)
SERVICE="http://redcap.nyumc.org/apps/redcap/api/"

# set the path to the REDCap R import script
REDCAP_IMPORT_R=/ifs/home/$(whoami)/projects/code/redcap_data_import_export.R
