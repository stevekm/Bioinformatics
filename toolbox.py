#!/usr/bin/env python


# a bunch of custom functions I written
# mostly for python 2.7, unless otherwise stated

# include this at the start of your Python script, with this file in the same dir:
# import toolbox as tb # my custom functions



def my_debugger(vars):
    # starts interactive Python terminal at location in script
    # call with tb.my_debugger(globals().copy()) anywhere in your script
    # or call my_debugger(locals().copy()) from anywhere within this package or another function
    import readline # optional, will allow Up/Down/History in the console
    import code
    # vars = globals().copy() # in python "global" variables are actually module-level
    vars.update(locals())
    shell = code.InteractiveConsole(vars)
    shell.interact()

def mkdir_p(path, return_path=False):
    # make a directory, and all parent dir's in the path
    import sys
    import os
    import errno

    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
    if return_path:
        return path

def initialize_file(string, output_file):
    # write string to file
    # !! THIS WILL OVERWRITE CONTENTS !!
    with open(output_file, "w") as myfile:
        myfile.write(string + '\n')

def append_string(string, output_file):
    # append string to file
    with open(output_file, "a") as myfile:
        myfile.write(string + '\n')

def subprocess_cmd(command):
    # run a terminal command with stdout piping enabled
    import subprocess as sp
    process = sp.Popen(command,stdout=sp.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    print proc_stdout


def get_files(dir_path, ends_with = '', trunc = False):
    # get the files in the dir that match the end pattern
    # trunc : return just the file dirname + basename (truncate)
    import sys
    import os
    file_list = []
    for subdir, dirs, files in os.walk(dir_path):
        for file in files:
            if file.endswith(ends_with):
                file_path = os.path.join(subdir,file)
                if (trunc):
                    file_dir = os.path.basename(os.path.dirname(file_path))
                    file_base = os.path.basename(file_path)
                    file_path = os.path.join(file_dir,file_base)
                file_list.append(file_path)
    return file_list


def download_file(my_URL, my_outfile = ''):
    # function to download a file from a URL
    # !! This will overwrite the output file
    # https://gist.github.com/hughdbrown/c145b8385a2afa6570e2

    import urllib2
    import urlparse
    import os

    URL_basename = os.path.basename(urlparse.urlsplit(my_URL).path)

    # if no output file specified, save to URL filename in current dir
    if my_outfile == '':
        my_outfile = URL_basename

    my_URL = urllib2.urlopen(my_URL)
    with open(my_outfile, 'wb') as output:
        while True:
            data = my_URL.read(4096) # download in chunks
            if data:
                output.write(data)
            else:
                break



def py_unzip(zip_file, outdir = "."):
    zip_ref = zipfile.ZipFile(zip_file, 'r')
    zip_ref.extractall(outdir)
    zip_ref.close()

def gz_unzip(gz_file, outdir = '', outfile = '', return_path = False):
    import gzip
    import os
    # extract a .gz file
    # !! This reads the entire file into memory !!

    # make sure the input file is a .gz file
    if not gz_file.lower().endswith('.gz'):
        print "ERROR: File is not a .gz file; ", gz_file
        return

    # read in the contents
    input_file = gzip.GzipFile(gz_file, 'rb')
    file_contents = input_file.read()
    input_file.close()

    # set the output path
    output_file_path = os.path.splitext(gz_file)[0]

    # if an outdir was passed, save the output there instead
    if outdir != '':
        output_file_path = os.path.join(outdir, os.path.basename(output_file_path))

    # if an output file was passed, use that instead
    if outfile != '':
        output_file_path = outfile

    # write the contents
    output_file = open(output_file_path, 'wb')
    output_file.write(file_contents)
    output_file.close()

    # return the path if requested
    if return_path:
        if os.path.exists(output_file_path):
            return output_file_path


def dict_from_tabular(inputfile, sep = ','):
    import csv
    lines_dict = {}
    reader = csv.reader(open(inputfile, 'r'), delimiter=sep)
    for key, value in reader:
        lines_dict[key] = value
    return lines_dict


def list_file_lines(file_path):
    # return the list of entries in a file, one per line
    # not blank lines, no trailing \n
    with open(file_path, 'r') as f:
        entries = [line.strip() for line in f if line.strip()]
    return entries


def split_df_col2rows(dataframe, split_col, split_char, new_colname):
    # # Splits a column into multiple rows
    # dataframe : pandas dataframe to be processed
    # split_col : chr string of the column name to be split
    # split_char : chr to split the col on
    # new_colname : new name for the
    # ~~~~~~~~~~~~~~~~ #
    import pandas as pd
    import numpy as np

    # make sure that the split_col is an 'object' type so we can split it
    if split_col in dataframe.select_dtypes([np.object]).columns:
        # save the split column as a separate object
        tmp_col = dataframe[split_col].str.split(split_char).apply(pd.Series, 1).stack()
        # drop the last index level
        tmp_col.index = tmp_col.index.droplevel(-1)
        # set the new col name
        tmp_col.name = new_colname
        # remove the original column from the df
        del dataframe[split_col]
        # join them into a new df
        dataframe = dataframe.join(tmp_col)
    else:
        print """
WARNING: Trying to split column {} in dataframe, where column is not dtype 'object'
Column dtype is: {}
Column will not be split but column name {} will be changed to: {}
        """.format(split_col, dataframe[split_col].dtype, split_col, new_colname)
        # just change the column name and keep moving
        dataframe.rename(columns={split_col: new_colname}, inplace=True)
    # rest the indexes
    dataframe = dataframe.reset_index(drop=True)
    return dataframe


def split_df_col2cols(dataframe, split_col, split_char, new_colnames, delete_old = False):
    # # Splits a column into multiple columns
    # dataframe : pandas dataframe to be processed
    # split_col : chr string of the column name to be split
    # split_char : chr to split the col on
    # new_colnames : list of new name for the columns
    # delete_old : logical True / False, remove original column?
    # ~~~~~~~~~~~~~~~~ #
    import pandas as pd
    import numpy as np
    # pl.my_debugger(globals().copy())
    # my_debugger(locals().copy())
    # save the split column as a separate object
    new_cols = dataframe[split_col].astype(np.object_).str.split(split_char).apply(pd.Series, 1)
    # if all values were NaN, no split occured, only one col exists still
    if len(new_cols.columns) < len(new_colnames):
        # create the missing cols, fill with NaN
        for i in range(len(new_cols.columns), len(new_colnames)):
            new_cols[new_colnames[i]] = np.nan
    # rename the cols
    new_cols.columns = new_colnames
    # remove the original column from the df
    if delete_old is True:
        del dataframe[split_col]
    # merge with df
    new_df = dataframe.join(new_cols)
    return new_df

def conjunction(*conditions):
    # apply multiple filtering conditions to a dataframe
    import numpy as np
    import functools
    return functools.reduce(np.logical_and, conditions)

def table_multi_filter(dataframe, filter_criteria):
    # filter a dataframe based on multiple criteria
    # 'filter_criteria' = {'include': {'column_name': ['value1', 'value2']}, ... }
    import pandas as pd
    # my_debugger(locals().copy())
    conditions_dfs = [] # empty list to hold conditional df's
    for key, value in filter_criteria['include'].items():
        for item in value:
            includes_df = dataframe[key] == item
            conditions_dfs.append(includes_df)

    for key, value in filter_criteria['exclude'].items():
        for item in value:
            excludes_df = dataframe[key] != item
            conditions_dfs.append(excludes_df)

    for key, value in filter_criteria['less_than'].items():
        less_thans_df = dataframe[key] < value
        conditions_dfs.append(less_thans_df)

    for key, value in filter_criteria['greater_than'].items():
        greater_thans_df = dataframe[key] > value
        conditions_dfs.append(greater_thans_df)

    for key, value in filter_criteria['less_or_null'].items():
        less_null_df = ( (dataframe[key] < value) | pd.isnull(dataframe[key]) )
        conditions_dfs.append(less_null_df)

    dataframe = dataframe[conjunction(*conditions_dfs)]
    return dataframe

def write_json(object, output_file):
    import json
    with open(output_file,"w") as f:
        json.dump(object, f, indent=4)

def load_json(input_file):
    import json
    with open(input_file,"r") as f:
        my_item = json.load(f)
    return my_item
