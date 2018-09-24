#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Divide a column in a table and output the value in a new column
"""
import os
import sys
import argparse

from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)
"""
https://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python
"""

def main(**kwargs):
    """
    """
    input_file = kwargs.pop('input_file', None)
    output_file = kwargs.pop('output_file', None)
    delim = kwargs.pop('delim', '\t')
    header = kwargs.pop('header', None)
    field_num = kwargs.pop('field_num')
    divisor = kwargs.pop('divisor')

    if input_file:
        fin = open(input_file)
    else:
        fin = sys.stdin

    if output_file:
        fout = open(output_file, "w")
    else:
        fout = sys.stdout

    if header:
        old_header = next(fin).strip()
        new_header = old_header + delim + header + '\n'
        fout.write(new_header)

    for line in fin:
        _field_num = int(field_num) - 1
        dividend = float(line.split(delim)[int(_field_num)]) / float(divisor)
        new_line = line.strip() + delim + str(dividend) + '\n'
        fout.write(new_line)

    fout.close()
    fin.close()

def parse():
    """
    Parses script args
    """
    parser = argparse.ArgumentParser(description='Append a column of text to a file')
    parser.add_argument("-i", default = None, dest = 'input_file', help="Input file")
    parser.add_argument("-o", default = None, dest = 'output_file', help="Output file")
    parser.add_argument("-d", default = '\t', dest = 'delim', help="Delimiter")
    parser.add_argument("--header", default = None, dest = 'header', help="Header for the new column")
    parser.add_argument("-f", "--field_num", required=True, dest = 'field_num', help="Field or column number in the table to be divide")
    parser.add_argument("--divisor", required=True, dest = 'divisor', help="Value to be divided by")
    args = parser.parse_args()

    main(**vars(args))

if __name__ == '__main__':
    parse()
