#!/usr/bin/env python
# python 2.7
'''
This script will recursively search all subdirectories for a given file, and
output all subdirectories that do not contain the given file

USAGE: find_missing_files.py ~/pipeline/align/results/ alignments.bam

ex:
/pipeline/align/results/sample1/alignments.bam
/pipeline/align/results/sample2
/pipeline/align/results/sample3/alignments.bam

output:
/pipeline/align/results/
/dir/to/search/align/sample2
'''
import sys
import os

search_this_dir = sys.argv[1]
missing_filename = sys.argv[2]

for dirpath, dirnames, filenames in os.walk(search_this_dir):
    if not missing_filename in filenames:
        print(dirpath)
