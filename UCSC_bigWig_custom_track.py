#!/usr/bin/python

import os
import glob
import pandas as pd
import csv

# ~~~~ GET SCRIPT ARGS ~~~~~~ #
# print 'Argument List:', str(sys.argv)
# print 'script name is:', sys.argv[0]

# ~~~~ COMMON ITEMS ~~~~~~ #
# location of the pipeline output
pipeline_align_dir="/ifs/home/kellys04/projects/SmithLab_ChIPSeq_RNASeq_2016_06_01/pipeline/align/results"

# location of this project output
analysis_outdir="/ifs/home/kellys04/projects/SmithLab_ChIPSeq_RNASeq_2016_06_01/project_notes/Combined_ChIP_RNA_Seq/analysis_output"

# external results dir for the bigwigs
external_bigwig_dir="/ifs/data/sequence/results/external/SmithLab/2016-06-08/Combined_RNA-Seq_ChIP-Seq/bigwig"

# http url to the external dir
external_URL="https://internet.edu/results/external/SmithLab/2016-06-08/Combined_RNA-Seq_ChIP-Seq/bigwig/"

# settings for the custom track views
bigwig_view_settings='visibility=full autoScale=off alwaysZero=on maxHeightPixels=50 graphType=bar viewLimits=0:1.0'


# ~~~~ GET FILES ~~~~~~ #
# find the bigwig files in the external dir
os.chdir(external_bigwig_dir)
# for file in glob.glob("*.bw"):
#     print(file)
bigWig_files=[file for file in glob.glob("*.bw")]


# ~~~~ CUSTOM TRACK FUNCITON ~~~~~~ #
def ucsc_cust_track_bigWig(track_name, track_url, view_settings):
	# this doesn't work in python 2.6.6: 
	# track_template='track type=bigWig name="{}" bigDataUrl={} {}'.format(track_name, track_url, view_settings) 
	track_template='track type=bigWig name="%s" bigDataUrl=%s%s %s' % (track_name, track_url, track_name, view_settings)
	
	return track_template


# ~~~~ MAKE TRACKS ~~~~~~ #

# for a column-wise per-sample set of tracks
# make an empty data frame to hold the sample outputs
track_table = pd.DataFrame()

for file in bigWig_files:
	print file
	
	# make the track
	tmp_track = ucsc_cust_track_bigWig(file,external_URL,bigwig_view_settings)
	
	# add the track to the dataframe in a new column
	track_table[file] = pd.Series([tmp_track])
	# in development, add more tracks for the sample.. 
	# index = name of track type

# print track_table to file
track_table.to_csv(analysis_outdir + '/UCSC_custom_track-wide.tsv',sep='\t', quoting=csv.QUOTE_NONE)

# for long table single column of tracks
track_table_long = pd.DataFrame()
track_table_long["col1"] = pd.Series([ucsc_cust_track_bigWig(file,external_URL,bigwig_view_settings) for file in bigWig_files])
# write to the file
track_table_long.sort(["col1"]).to_csv(analysis_outdir + '/UCSC_custom_track-long.tsv',sep='\t',index=False,header=False, quoting=csv.QUOTE_NONE)

# TO-DO: add more types of tracks! see the bash version of this for functionality to port over to Python
