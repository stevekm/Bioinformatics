#!/usr/bin/env python
# python 2.7

'''
USAGE: dump_xlsx.py /path/to/my_file.xlsx

OUTPUT: /path/to/my_file.sheet_1.tsv

This script will dump every sheet in an XLSX Excel file to a TSV
'''

import pandas as pd
import numpy as np
import sys
import os


file_name = sys.argv[1]
file_base = os.path.splitext(os.path.basename(file_name))[0]
file_dir = os.path.dirname(file_name)


# read excel file
xls_file = pd.ExcelFile(file_name)

# load each excel sheet into a dict entry
xls_dict = {sheet_name: xls_file.parse(sheet_name) for sheet_name in xls_file.sheet_names}

# read each sheet into a Pandas dataframe then export as TSV
count = 1
for sheet_name, sheet_data, in xls_dict.iteritems():
    sheet_df = xls_dict[sheet_name]
    out_file_name = '.'.join([file_base, sheet_name, "tsv"]).replace(" ", "_") # use this if you want to preserve the sheet names
    # out_file_name = '.'.join([file_base, "sheet_" + str(count), "tsv"]) # use this if you don't want to preserve sheet names
    out_file_path = os.path.join(file_dir, out_file_name)
    sheet_df.to_csv(out_file_path,sep ='\t', index = False, encoding = 'utf-8')
    count += 1
