#!/usr/bin/env python
# python 2.7

'''
This script will attempt to convert the massive methylation data set
which is saved in a 100000000000000MB XLSX file
to a text format so it can be more easily parsed
'''

import pandas as pd
import numpy as np
import sys
import os
import pickle

def save_pydata(data, outfile):
    # save py data in pickle format
    with open(outfile, 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
        print 'Object saved to file:\n', outfile


def load_pydata(infile):
    # open py pickle data
    with open(infile, 'rb') as f:
        # The protocol version used is detected automatically, so we do not
        # have to specify it.
        data = pickle.load(f)
    return data


# file_name = sys.argv[1]
file_name = "methylationData.ratio.norm.178_L2_with_annotation_TB.xlsx"

# read excel file
xls_file = pd.ExcelFile(file_name)

# load each excel sheet into a dict entry
xls_dict = {sheet_name: xls_file.parse(sheet_name) for sheet_name in xls_file.sheet_names}
save_pydata(xls_dict, "methyl_data_pandas_xls_dict.pickle")


# get the first sheet dataframe
# xls_df = xls_dict.itervalues().next()
# print xls_df.head()

# >>> print xls_dict.keys()
# [u'Sheet1', u'methylationData.ratio.norm.178_']

sheet1_df = xls_dict['Sheet1']
save_pydata(sheet1_df, "methyl_data_pandas_sheet1_df.pickle")
sheet1_df.to_csv("methyldata_sheet1.tsv",sep ='\t', index = False)

methyldata_df = xls_dict['methylationData.ratio.norm.178_']
save_pydata(methyldata_df, "methyl_data_pandas_methyldata_df.pickle")
methyldata_df.to_csv("methyl_data_ratio_norm_178.tsv",sep ='\t', index = False)

