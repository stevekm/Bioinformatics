#!/usr/bin/env python
# python 2.7
'''
This script will output directories which are missing any of the given files

USAGE: find_missing_files.py ~/pipeline/align/results/ alignments.bam track.bw
'''
import sys
import os

# get a copy of the script args
copy_args = list(sys.argv)
# remove the scriptname
copy_args.pop(0)

# search location = first script arg
search_this_dir = copy_args.pop(0)
# all other script args = filenames to search for
missing_filenames = copy_args


for dirpath, dirnames, filenames in os.walk(search_this_dir):
    # check if the dir is missing any of the desired files
    if not set(missing_filenames).issubset(filenames):
        print(dirpath)
