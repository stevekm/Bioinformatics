#/usr/bin/env python

'''
This script will run gtools in order to generate coverage matrices for a large number of bedgraphs with corresponding region BED files
'''


import os
import sys
import csv
import toolbox as tb # my custom functions


class projectLocs:
    '''
    An object class to hold locations of places in the project
    '''
    def __init__(self):
        self.project_dir = "/ifs/home/kellys04/projects/Smithlab_ChIpSeq_2016-03-10/project_notes/methylation_profiles"
        self.bedgraphs_dir = os.path.join(self.project_dir, "bedgraphs")
        self.regions_dir = os.path.join(self.project_dir, "sample_gene_TSS_5Kbp_regions")
        self.matrix_dir = os.path.join(self.project_dir, "matrices")
        self.matrix_logdir = tb.mkdir_p(os.path.join(self.matrix_dir, "logs"), return_path=True)
        self.samplesheet = os.path.join(self.project_dir, "microarray_methlyation_samplesheet_3.tsv")


def file_match(dir, start_pattern = '', end_pattern = '', contains_pattern = ''):
    '''
    Find a file in a dir which matches the supplied patterns
    NOTE: Doesn't search recursively!
    '''
    file_match = []
    for file in os.listdir(dir):
        if file.startswith(start_pattern) and file.endswith(end_pattern) and contains_pattern in file:
            file_match.append(os.path.join(dir, file))
    return file_match


def samplesheet2dict(samplesheet_file, sample_colname, sep = '\t'):
    '''
    Create a nested dict for each sample in a standard TSV samplesheet with headers
    NOTE: If duplicate column headers exist, latter will overwrite former
    '''
    import csv
    sample_files_dict = {}
    with open(samplesheet_file) as csvfile:
        reader = csv.DictReader(csvfile, delimiter = sep)
        for row_dict in reader:
            sample_files_dict[row_dict[sample_colname]] = row_dict
    return sample_files_dict

def print_samplesheet_dict(samplesheet_dict):
    '''
    Pretty printing for a samplesheet dict so its more readable on the console
    '''
    for key, value in samplesheet_dict.iteritems():
        print key
        for subkey, subvalue in value.iteritems():
            print subkey,':', subvalue
        print "\n"


# sample sheet looks like this:
# Sample_Name	Methylation_name	Methylation_name_D	Methylation_name_R	Microarray_name_D	Microarray_name_R	Microarray_name_LR	D_Status	R_Status	Microarray_name_D	Microarray_name_R	Microarray_name_LR
# AGK	PAPAGK	PAPAGK.Diagnosis	PAPAGK.Relapse	AGK-D	AGK-R	AGK-LR	D	R	AGK-Exp-D	AGK-Exp-R	AGK-Exp-LR


# create object to hold the project locations
proj_locs = projectLocs()

# create the samplesheet dict
samplesheet_dict = samplesheet2dict(samplesheet_file = proj_locs.samplesheet, sample_colname = 'Sample_Name')

# print it to the console
# print_samplesheet_dict(samplesheet_dict)

# samplesheet IDs whose values match filenames
regions_r_ID_pattern = 'Microarray_name_R' # SPN-R
regions_d_ID_pattern = 'Microarray_name_D'
methyl_r_ID_pattern = 'Methylation_name_R' # PAPSPN.Relapse
methyl_d_ID_pattern = 'Methylation_name_D'


# filename patterns
regions_top_pattern = 'Top_expressed_genes.bed'
regions_bottom_pattern = 'Bottom_expressed_genes.bed'
bedgraph_pattern = '.bedgraph'

# make lists of the items from above to iterate over
region_IDs = [regions_r_ID_pattern, regions_d_ID_pattern]
region_expressions = [regions_top_pattern, regions_bottom_pattern]
methyl_IDs = [methyl_r_ID_pattern, methyl_d_ID_pattern]


def qsub_gtools_matrix(samplesheet_dict, proj_locs, region_IDs, region_expressions, methyl_IDs, bedgraph_pattern):
    '''
    Submit a qsub job to run every gtools matrix on the combinations of bedgraphs and region files
    '''
    # parameters for qsub job
    job_threads = "1"
    job_mem = "4G"
    job_options = "-j y" # merge stderr and stdout # job_options="-l mem_free=$job_mem -l h_vmem=$job_mem -l mem_token=$job_mem"
    for sampleID, items in samplesheet_dict.iteritems():
        print(sampleID)
        for region_ID in region_IDs:
            print(samplesheet_dict[sampleID][region_ID])
            for region_expression in region_expressions:
                print(region_expression)
                for methyl_ID in methyl_IDs:
                    print(samplesheet_dict[sampleID][methyl_ID])
                    # find the BED file for the combination of Sample + genes expression
                    regions_file = file_match(proj_locs.regions_dir, start_pattern = samplesheet_dict[sampleID][region_ID], end_pattern = region_expression)[0]
                    # find the bedgraph file with the values from the methylation analysis
                    bedgraph_file = file_match(proj_locs.bedgraphs_dir, start_pattern = samplesheet_dict[sampleID][methyl_ID], end_pattern = bedgraph_pattern)[0]
                    # set up the output file naming scheme
                    output_file_base = '{}_{}_{}'.format(
                    samplesheet_dict[sampleID][region_ID],
                    region_expression,
                    samplesheet_dict[sampleID][methyl_ID]
                    )
                    output_file_basename = '{}.matrix'.format(output_file_base)
                    #
                    output_file = os.path.join(proj_locs.matrix_dir, output_file_basename)
                    # the command to run gtools
                    gtools_command = '''
set -x
head "{}"
head "{}"
gtools-threaded matrix -v -i -o {} -i -nbins 51 --overlap-op value -rpkm -profile sum {} {}
head "{}"
                    '''.format(bedgraph_file,
                    regions_file,
                    output_file,
                    bedgraph_file,
                    regions_file,
                    output_file)
                    # the command to submit gtools to the cluster; REQUIRES BASH !
                    qsub_command = '''
mkdir -p "{}" # make sure the log dir exists
qsub -wd "{}" -o :{}/ -e :{}/ -pe threaded {} -N "{}" {} <<E0F
{}
E0F
                    '''.format(
                    proj_locs.matrix_logdir,
                    proj_locs.project_dir,
                    proj_locs.matrix_logdir, proj_locs.matrix_logdir,
                    job_threads,
                    output_file_base,
                    job_options,
                    gtools_command)
                    #
                    # print(qsub_command)
                    # tb.my_debugger(globals().copy())
                    tb.subprocess_cmd(qsub_command)

# run the gtools qsub functions
qsub_gtools_matrix(samplesheet_dict, proj_locs, region_IDs, region_expressions, methyl_IDs, bedgraph_pattern)

sys.exit()
