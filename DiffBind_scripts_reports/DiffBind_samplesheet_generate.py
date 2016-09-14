#!/usr/bin/env python

# Generate a sample sheet for DiffBind analysis using the pipeline results
# using Python 2.7

import sys
import os
import pandas as pd
import subprocess
import errno

# use terminal 'find' to return file path
def get_file_path(search_dir='.',name='',path='', more=''):
    print 'Searching for file:\t', name, '\nIn location:\t', search_dir + path
    filepath = [line for line in subprocess.check_output("find {} -name '{}' -path '{}'{}".format(search_dir,name,path,more), shell=True).splitlines()]
    return filepath[0] # return first result

def mkdir_p(path):
    # if not os.path.exists(dir_path):
    # os.makedirs(dir_path)
    # distutils.dir_util.mkpath(name[, mode=0777, verbose=0, dry_run=0])
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

project_sample_sheet_file = "/ifs/data/smithlab/SmithLab-ChIPSeq_2016-06-06/inputs/sample-sheet.tsv"
# sample    control group   fastq-r1    fastq-r2    genome  fragmentation-size  cellline    treatment   chip    cellline_treatment  
align_results_dir="/ifs/data/smithlab/SmithLab-ChIPSeq_2016-06-06/pipeline/align/results"
peaks_results_dir="/ifs/data/smithlab/SmithLab-ChIPSeq_2016-06-06/pipeline/peaks/results" 
peak_caller="macs"
out_dir = "/ifs/home/kellys04/projects/SmithLab-ChIPSeq_2016-06-06/analysis_dir/project_notes/diffbind_analysis"
diffbind_outdir = out_dir + '/output'
diffbind_analysis_script = '/ifs/data/smithlab/SmithLab-ChIPSeq_2016-06-06/project_notes/diffbind_analysis/code/chipseq-diffbind.R'
diffbind_samplesheet_output = out_dir + '/diffbind_samplesheet2.csv'

mkdir_p(diffbind_outdir)

# read in the project sample sheet
sample_sheet = pd.read_table(project_sample_sheet_file,sep='\t',header=0)

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
    celllineID = sample_sheet.loc[i,'cellline']
    ChIPID = sample_sheet.loc[i,'chip']
    treatmentID = sample_sheet.loc[i,'treatment']
    if not pd.isnull(controlID):
        sample_bam = get_file_path(align_results_dir,'*.bam','*{}*'.format(sampleID))
        control_bam = get_file_path(align_results_dir,'*.bam','*{}*'.format(controlID))
        sample_bed = get_file_path(peaks_results_dir,'macs_peaks.xls','*{}*'.format(sampleID), " -path '*peaks.by_sample.macs_narrow/align.by_sample.bowtie2*'")
        dict1.update({'SampleID':sampleID, 'Tissue':celllineID, 'Factor':ChIPID, 'Condition':conditionID, 'Treatment':treatmentID, 'Replicate':'', 'bamReads':sample_bam, 'bamControl':control_bam, 'Peaks':sample_bed, 'PeakCaller':peak_caller})
        rows_list.append(dict1)

# make data frame
df = pd.DataFrame(rows_list, columns=diffbind_sheet_columns)

print '\nDiffBind Sample sheet:\n'
print df

# save to file
print '\nSample sheet saved to:\n', diffbind_samplesheet_output, '\n'
df.to_csv(diffbind_samplesheet_output, index = False)


