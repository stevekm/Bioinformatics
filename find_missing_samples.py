#!/usr/bin/env python
# python 2.7
'''
USAGE: find_missing_samples.py /path/to/sample-sheet.tsv /path/to/results_dir

EXAMPLE: find_missing_samples.py ~/projects/inputs/sample-sheet.tsv ~/projects/pipeline/align/results/

This script will check the output directory of the analysis pipeline to make sure
all the samples listed in the sample sheet are present in the output.

For every sample in the samplesheet, there should be a corresponding directory
in the results dir which has no subdirectories, and has the same name as a sample ID

This script will find the sampleIDs that do not have a corresponding results subdir, and vice versa
'''
import sys
import os
import csv

samplesheet_file = sys.argv[1] # should be tab separated!
results_dir = sys.argv[2]

# get sample IDs from the first column of the samplsheet
sampleIDs = []
with open(samplesheet_file) as ss:
    header_line = next(ss)
    csvreader = csv.reader(ss, delimiter='\t')
    for line in csvreader:
        if len(line) > 0: # skip empty lines
            sampleIDs.append(line[0])


# find the lowest dirs in the results directory; dirnames = sample IDs
dir_sampleIDs = []
for dirpath, dirnames, filenames in os.walk(results_dir):
    if not dirnames:
        dir_basename = os.path.basename(dirpath)
        # skip the .db dir created by the pipeline
        if dir_basename != '.db':
            dir_sampleIDs.append(dir_basename)

# sampleIDs with no corresponding results dir
missing_items = [x for x in sampleIDs if x not in dir_sampleIDs]
if len(missing_items) > 0:
    print "The following sample IDs were not found in the results dir:"
    for item in missing_items:
        print item

# results dirs with no corresponding sample ID
mystery_dirs = [x for x in dir_sampleIDs if x not in sampleIDs]
if len(mystery_dirs) > 0:
    print "The following results dirs did not match samplesheet IDs:"
    for item in mystery_dirs:
        print item
