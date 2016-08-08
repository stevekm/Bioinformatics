#!/usr/bin/env python

# This workflow will search for all VCF files in the given directory, 
# run ANNOVAR on them, 
# and then process and merge the annotations to a single file
# with labels per original file & Sample ID per variant
# # developed with python 2.7

# qual > 100 ; VCF
# coverage > 250 ; VCF
# VF > 5 ; VCF
# Exonic ; ANNOVAR
# Nonsynonmous ; ANNOVAR
# MAF < 0.01 ; ANNOVAR
# SB < 0.8 ; VCF
# annotate, filter, ID, merge all


# ~~~~ LOAD PACKAGES ~~~~~~ #
import sys
import os
import pandas as pd
import numpy as np
import fnmatch
import re
import subprocess as sp
import csv


# ~~~~ GET SCRIPT ARGS ~~~~~~ #
# python test.py arg1 arg2 arg3
# print 'Number of arguments:', len(sys.argv), 'arguments.'
# print 'Argument List:', str(sys.argv)
# print 'script name is:', sys.argv[0]

# ~~~~ COMMON LOCATIONS & ITEMS ~~~~~~ #
project_dir = "/ifs/home/kellys04/projects/variant_reporting"
input_dir = (project_dir + "/input")
output_dir = (project_dir + "/output")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)


# ~~~~ CUSTOM FUNCTIONS ~~~~~~ #
def subprocess_cmd(command):
    # proc = sp.Popen(['ls','-l'])
    process = sp.Popen(command,stdout=sp.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    print proc_stdout

def convert_2_annovar(vcf_input_file, av_output_file):
    if not os.path.exists(av_output_file) and not os.path.isfile(av_output_file):
        subprocess_cmd('convert2annovar.pl -format vcf4old {} -includeinfo > {}'.format(vcf_input_file,av_output_file))


def annovar_table(avinput, annovar_output, annovar_db_dir="/ifs/home/kellys04/software/annovar/db", build_version="hg19", annovar_protocol="-protocol refGene,cosmic68,clinvar_20150629,1000g2015aug_all -operation g,f,f,f"
):
    # avinput : file path to avinput
    # annovar_output : file path to annovar output
    #  --otherinfo # http://annovar.openbioinformatics.org/en/latest/user-guide/startup/
    file_suffix = '.' + build_version + '_multianno.txt'
    output_file = annovar_output + file_suffix
    if not os.path.exists(output_file) and not os.path.isfile(output_file):
        subprocess_cmd('table_annovar.pl {} {} -buildver {} -out {} -remove {} -nastring .'.format(avinput, annovar_db_dir, build_version, annovar_output, annovar_protocol))
    # file.endswith('multianno.txt') # TSVC_variants_IonXpress_015_filtered.hg19_multianno.txt
    


def vcf_header_skip(file_path):
    skip_rows = 0
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('##'):
                skip_rows += 1
            else:
                break
    return skip_rows

def vcf_qual_filter(vcf_table, qual_score = 100):
    # QUAL > 100 # df.loc[df['column_name'] == some_value]
    filtered_table = vcf_table.loc[raw_vcf_table['QUAL'] > qual_score]
    return filtered_table


def av_table_filter(av_table, gene_func = "exonic", exonic_func = "nonsynonymous SNV", maf_cutoff = 0.01):
    # Func.refGene : exonic
    # ExonicFunc.refGene : nonsynonymous SNV
    # 1000g2015aug_all : < 0.01
    filtered_table = av_table.loc[ (av_table['ExonicFunc.refGene'] == exonic_func) & (av_table['1000g2015aug_all'] < maf_cutoff) & (av_table['Func.refGene'] == gene_func) ]
    return filtered_table




# ~~~~ ANNOTATION PIPELINE ~~~~~~ #
# ~~~~ FIND VCF FILES ~~~~~~ #
# empty list to hold filepaths
raw_vcf_paths = []
filtered_vcf_paths = []
annovar_output_paths = []

# find the VCF files in the input dir
for subdir, dirs, files in os.walk(input_dir):
    for file in files:
    	if file.endswith('.vcf'):
			raw_vcf_paths.append((os.path.join(subdir,file)))

raw_vcf_paths.sort()



# ~~~~ FILTER VCF QUALITY ~~~~~~ #
filtered_vcf_dir = output_dir + '/filtered_vcf'
if not os.path.exists(filtered_vcf_dir):
    os.makedirs(filtered_vcf_dir)

# read in the file with the row skips
for file in raw_vcf_paths:
    filtered_vcf_file = filtered_vcf_dir + '/' + os.path.basename(file).replace('.vcf', '_filtered.vcf')
    raw_vcf_table = pd.read_table(raw_vcf_paths[1],sep='\t',skiprows=vcf_header_skip(file),header=0)
    filtered_vcf_table = vcf_qual_filter(raw_vcf_table)
    filtered_vcf_table.to_csv(filtered_vcf_file, sep='\t', index=False)
    filtered_vcf_paths.append(filtered_vcf_file)



# ~~~~ ANNOTATE VCF FILES ~~~~~~ #
file_paths = filtered_vcf_paths # tmp for testing
for file in file_paths:
    file_basename = os.path.basename(file)
    avinput = file.replace('.vcf', '.avinput')
    annovar_output = output_dir + '/' + file_basename.replace('.vcf', '') 
    convert_2_annovar(file,avinput)    
    annovar_table(avinput, annovar_output)




# ~~~~ FILTER ANNOVAR OUTPUT ~~~~~~ #
# find the files
for subdir, dirs, files in os.walk(output_dir):
    for file in files:
        if file.endswith('multianno.txt'):
            annovar_output_paths.append((os.path.join(subdir,file)))

# filter files

for file in annovar_output_paths:
    # print file
    av_table = pd.read_table(file)
    print av_table_filter(av_table)



sys.exit()



