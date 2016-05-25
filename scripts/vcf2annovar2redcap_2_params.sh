#!/bin/bash

# params file for the vcf2annovar2redcap_2 script
# location of the annovar convert program
ANNOVAR_CONVERT=/path/to/software/annovar/convert2annovar.pl

# table command
ANNOVAR_TABlE=/path/to/software/annovar/table_annovar.pl

# location of the ANNOVAR databases; I changed the name of the dir btw
DB_DIR=/path/to/software/annovar/db

# Convert the VCF to TSV, using libvcf <https://github.com/ekg/vcflib>
VCF_to_TSV=~/software/vcflib/bin/vcf2tsv

# location of the ANNOVAR conversion script to produce REDCap output
ANNOVAR2REDCAP=/path/to/code/annovar2redcap.R

# ANNOVAR build version, e.g. hg19
BUILD_VERSION=hg19
