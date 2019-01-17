#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Reformats the .tsv vcf table output to recalculate and standardize values for downstream processing.

INPUT: .TSV formatted output from GATK VariantsToTable
OUTPUT: .TSV formatted table with recalculated variant allele frequency in the FREQ column
USAGE: reformat-vcf-table.py -c GATKHC -s "sampleID" -i "sample_tsv" -o "sampleID.recalc.tsv"
"""
import csv
import sys
import argparse

from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)
"""
https://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python
"""

def GATKHC(fin, fout, sampleID):
    """
    Recalculates the variant allele frequency for GATK Haplotype Caller output
    assumes VCFv4.2
    Outputs extra columns using variant values

    Parameters
    ----------
    fin: file connection
        connection to the input file
    fout: file connection
        connection to the output file
    sampleID: str
        identifier for the sample in the input file connection

    """
    # column names for the AD and DP fields in the table
    AD_key = "{0}.AD".format(sampleID)
    DP_key = "{0}.DP".format(sampleID)

    reader = csv.DictReader(fin, delimiter = '\t')
    fieldnames = reader.fieldnames
    fieldnames.append('FREQ')
    fieldnames.append('AD.REF')
    fieldnames.append('AD.ALT')
    fieldnames.append('AF.ALT')
    fieldnames.append('AF.REF')
    fieldnames.append('DP')
    fieldnames.append('AF')
    fieldnames_out = [ item for item in fieldnames ]
    fieldnames_out.remove(AD_key)
    fieldnames_out.remove(DP_key)
    writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = fieldnames_out)
    writer.writeheader()
    for row in reader:
        ref_AD = float(row[AD_key].split(',')[0])
        alt_AD = float(row[AD_key].split(',')[1])
        depth = float(row[DP_key])
        ref_AF = ref_AD / depth
        alt_AF = alt_AD / depth
        row['FREQ'] = alt_AF
        row['AF'] = alt_AF
        row['AD.REF'] = ref_AD
        row['AD.ALT'] = alt_AD
        row['AF.REF'] = ref_AF
        row['AF.ALT'] = alt_AF
        row['DP'] = depth
        drop_vals = []
        drop_vals.append(row.pop(AD_key))
        drop_vals.append(row.pop(DP_key))
        writer.writerow(row)

def LoFreq(fin, fout):
    """
    LoFreq does not need recalulating; output a new column in the file with 'FREQ' for consistency

    Parameters
    ----------
    fin: file connection
        connection to the input file
    fout: file connection
        connection to the output file

    """
    # allele frequency column
    AF_key = "AF"
    reader = csv.DictReader(fin, delimiter = '\t')
    fieldnames = reader.fieldnames
    fieldnames.append('FREQ')
    writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = fieldnames)
    writer.writeheader()
    for row in reader:
        AF_val = row[AF_key]
        row['FREQ'] = AF_val
        writer.writerow(row)

def MuTect2(fin, fout):
    """
    Outputs the TUMOR.AF as a new column called 'FREQ', and outputs the TLOD value as QUAL
    Adds extra allelic depth columns

    Parameters
    ----------
    fin: file connection
        connection to the input file
    fout: file connection
        connection to the output file

    """
    reader = csv.DictReader(fin, delimiter = '\t')
    # get old headers
    fieldnames = reader.fieldnames
    # append new headers for the columns to be created
    fieldnames.append('FREQ')
    fieldnames.append('QUAL')
    fieldnames.append('TUMOR.AD.REF')
    fieldnames.append('TUMOR.AD.ALT')
    fieldnames.append('TUMOR.AD.TOTAL')
    fieldnames.append('NORMAL.AD.REF')
    fieldnames.append('NORMAL.AD.ALT')
    fieldnames.append('NORMAL.AD.TOTAL')
    writer = csv.DictWriter(fout, delimiter = '\t', fieldnames = fieldnames)
    writer.writeheader()
    for row in reader:
        # set the FREQ to the tumor AF
        tumor_AF_value = row["TUMOR.AF"]
        row['FREQ'] = row['TUMOR.AF']
        # change QUAL to the TLOD
        row['QUAL'] = row['TLOD']
        # split the tumor AD values for ref and alt
        row['TUMOR.AD.REF'] = row['TUMOR.AD'].split(',')[0]
        row['TUMOR.AD.ALT'] = row['TUMOR.AD'].split(',')[1]
        row['NORMAL.AD.REF'] = row['NORMAL.AD'].split(',')[0]
        row['NORMAL.AD.ALT'] = row['NORMAL.AD'].split(',')[1]
        # add up the total allelic depths
        row['TUMOR.AD.TOTAL'] = int(row['TUMOR.AD.REF']) + int(row['TUMOR.AD.ALT'])
        row['NORMAL.AD.TOTAL'] = int(row['NORMAL.AD.REF']) + int(row['NORMAL.AD.ALT'])
        writer.writerow(row)

def main(**kwargs):
    """
    Main control function for the script
    """
    input_file = kwargs.pop('input_file', None)
    output_file = kwargs.pop('output_file', None)
    caller = kwargs.pop('caller')
    sampleID = kwargs.pop('sampleID')

    if input_file:
        fin = open(input_file)
    else:
        fin = sys.stdin

    if output_file:
        fout = open(output_file, "w")
    else:
        fout = sys.stdout

    if caller == "GATKHC":
        GATKHC(fin, fout, sampleID)
        fout.close()
        fin.close()
    elif caller == "LoFreq":
        LoFreq(fin, fout)
        fout.close()
        fin.close()
    elif caller == "MuTect2":
        MuTect2(fin, fout)
        fout.close()
        fin.close()
    else:
        print("ERROR: caller not recognized: {0}".format(caller))
        sys.exit(1)




def parse():
    """
    Parses script args
    """
    parser = argparse.ArgumentParser(description='Append a column of text to a file')
    parser.add_argument("-i", default = None, dest = 'input_file', help="Input file")
    parser.add_argument("-o", default = None, dest = 'output_file', help="Output file")

    parser.add_argument("-c", "--caller", dest = 'caller', help="Variant caller used", required=True)
    parser.add_argument("-s", "--sampleID", dest = 'sampleID', help="Sample ID", required=True)
    args = parser.parse_args()

    main(**vars(args))



if __name__ == '__main__':
    parse()
