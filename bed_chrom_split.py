#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Split a .bed file into a single file for each chromosome
'''
import os

# ~~~~~ FUNCTIONS ~~~~~ #
def get_bed_chroms(bed_file):
    '''
    Return a list of all chromosomes in a .bed file
    '''
    chroms = []
    with open(bed_file) as f:
        for line in f:
            parts = line.split()
            # make sure there's an entry
            if len(parts[0]) > 0:
                chroms.append(parts[0])
    # get unique entries
    chroms = list(set(chroms))
    return(chroms)


def make_bed_splitchrom_filenames(bed_file, chroms = None, output_dir = None):
    '''
    Make a dict of filenames for a per-chromosome split .bed file based on
    the chroms present

    ex:
    {'chr15': 'targets_chr15.bed', 'chr14': 'targets_chr14.bed'}
    '''
    if not chroms:
        chroms = get_bed_chroms(bed_file)
    if output_dir:
        file_dirpath = output_dir
    else:
        file_dirpath = os.path.dirname(bed_file)

    file_dirpath = os.path.realpath(file_dirpath)

    file_base, file_ext = os.path.splitext(os.path.basename(bed_file))
    chrom_filenames = {}
    for chrom in chroms:
        filename = '{0}_{1}{2}'.format(file_base, chrom, file_ext)
        file_path = os.path.join(file_dirpath, filename)
        chrom_filenames.update({chrom: file_path})
    return(chrom_filenames)


def append_line(output_file, line):
    '''
    Append a line of text to a file
    '''
    with open(output_file, 'a') as f:
        f.write(line)

def split_bed_by_chrom(bed_file, **kwargs):
    '''
    Split a .bed file into sub-files for each chromosome present

    split_bed_by_chrom(bed_file = bed_file, output_dir = output_dir)
    '''
    # get the full path to the file
    bed_file = os.path.realpath(bed_file)
    chroms = get_bed_chroms(bed_file)
    chrom_filenames = make_bed_splitchrom_filenames(bed_file = bed_file, chroms = chroms, **kwargs)
    # print(chrom_filenames)
    with open(bed_file) as f:
        for line in f:
            parts = line.split()
            if chrom_filenames.get(parts[0], None):
                chrom = parts[0]
                output_file = chrom_filenames[chrom]
                append_line(output_file = output_file, line = line)



# ~~~~~ RUN ~~~~~ #
bed_file = 'targets.bed'
output_dir = 'test'
split_bed_by_chrom(bed_file = bed_file, output_dir = output_dir)
