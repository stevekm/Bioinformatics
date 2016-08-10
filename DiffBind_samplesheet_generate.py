#!/usr/bin/env python

# Generate a sample sheet for DiffBind analysis using the pipeline results
# using Python 2.7

import os
import pandas as pd
import subprocess

# use terminal 'find' to return file path
def get_file_path(search_dir='.',name='',path='', more=''):
    print 'Searching for file:\t', name, '\nIn location:\t', search_dir + path
    filepath = [line for line in subprocess.check_output("find {} -name '{}' -path '{}'{}".format(search_dir,name,path,more), shell=True).splitlines()]
    return filepath[0] # return first result

sample_sheet_file = "/ifs/data/smithlab/SmithLab-ChIPSeq_2016-06-06/inputs/sample-sheet.tsv"
align_results_dir="/ifs/data/smithlab/SmithLab-ChIPSeq_2016-06-06/pipeline/align/results"
peaks_results_dir="/ifs/data/smithlab/SmithLab-ChIPSeq_2016-06-06/pipeline/peaks/results" 
peak_caller="macs"
out_dir = "/ifs/home/kellys04/projects/SmithLab-ChIPSeq_2016-06-06/analysis_dir/project_notes/diffbind_analysis"

# read in the project sample sheet
sample_sheet = pd.read_table(sample_sheet_file,sep='\t',header=0)

# columns to use the in DiffBind sample sheet
diffbind_sheet_columns = ['SampleID', 'Tissue', 'Factor', 'Condition', 'Treatment', 'Replicate', 'bamReads', 'bamControl', 'Peaks', 'PeakCaller']

# list of entries for the final DiffBind sample sheet dataframe
rows_list = []

# iterate over the pipeline sample sheet rows
print '\nGenerating DiffBind sample sheet\n'
for i in sample_sheet.index:
    dict1 = {}
    sampleID = sample_sheet.loc[i,'sample']
    controlID = sample_sheet.loc[i,'control']
    conditionID = sample_sheet.loc[i,'cellline_treatment']
    if not pd.isnull(controlID):
        sample_bam = get_file_path(align_results_dir,'*.bam','*{}*'.format(sampleID))
        control_bam = get_file_path(align_results_dir,'*.bam','*{}*'.format(controlID))
        sample_bed = get_file_path(peaks_results_dir,'macs_peaks.xls','*{}*'.format(sampleID), " -path '*peaks.by_sample.macs_narrow/align.by_sample.bowtie2*'")
        dict1.update({'SampleID':sampleID, 'Tissue':'?', 'Factor':'?', 'Condition':conditionID, 'Treatment':'?', 'Replicate':'', 'bamReads':sample_bam, 'bamControl':control_bam, 'Peaks':sample_bed, 'PeakCaller':peak_caller})
        rows_list.append(dict1)

# make data frame
df = pd.DataFrame(rows_list, columns=diffbind_sheet_columns)

print '\nDiffBind Sample sheet:\n'
print df

# save to file
print '\nSample sheet saved to:\n', out_dir + '/diffbind_samplesheet.csv\n'
df.to_csv(out_dir + '/diffbind_samplesheet.csv', index = False)


