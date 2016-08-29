#!/usr/bin/env python
# python 2.7

from Bio import SeqIO
import re

def search_pattern(pattern, file_path):
    num_match = 0
    line_number = 0
    with open(file_path, 'r') as f:
            for line in f:
                line_number += 1
                if pattern in line:
                    num_match += 1
                    print line_number, line
    print "Number of matches found:\t", num_match
