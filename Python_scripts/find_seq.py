#!/usr/bin/env python
# python 2.7

# find an exact match of a nucleotide sequence within a FASTA file
# USAGE: find_seq.py <pattern> /path/to/file.fasta
# example: find_seq.py TAACCACTCA /local/data/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa > test_loc.txt

from Bio import SeqIO
import re
import sys


# pattern = "TAACCACTCA"
pattern = sys.argv[1]

# file_path = "/local/data/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
file_path = sys.argv[2]

def search_fasta(pattern, file_path):
    for record in SeqIO.parse(open(file_path, "rU"), "fasta"):
        chrom = record.id
        for match in re.finditer(pattern, str(record.seq)):
            start_pos = match.start() + 1
            end_pos = match.end() + 1
            print chrom, '\t', start_pos, '\t', end_pos

search_fasta(pattern, file_path)

