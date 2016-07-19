#!/usr/bin/python

import os
import glob
import pandas as pd
import csv

# https://genome.ucsc.edu/goldenpath/help/customTrack.html

# ~~~~ GET SCRIPT ARGS ~~~~~~ #
# print 'Argument List:', str(sys.argv)
# print 'script name is:', sys.argv[0]

# ~~~~ COMMON ITEMS ~~~~~~ #
# external results dir for the bigwigs
chipseq_external_bigwig_dir="/ifs/data/sequence/results/external/smithlab/bigWigs/ChIP-Seq"
rnaseq_external_bigwig_dir="/ifs/data/sequence/results/external/smithlab/bigWigs/RNA-Seq"

# http url to the external dir
chipseq_external_URL="https://genome.med.nyu.edu/results/external/smithlab/bigWigs/ChIP-Seq/"
rnaseq_external_URL="https://genome.med.nyu.edu/results/external/smithlab/bigWigs/RNA-Seq/"

# settings for the custom track views
chipseq_bigwig_view_settings='visibility=full autoScale=off alwaysZero=on maxHeightPixels=50 graphType=bar viewLimits=0:1.5 color=0,0,255'
rnaseq_bigwig_view_settings='visibility=full autoScale=off alwaysZero=on maxHeightPixels=50 graphType=bar viewLimits=0:3 color=255,0,0'
# chr16:8,968,103-9,075,189
# USP7

# location for the output from this script
analysis_outdir="/ifs/home/kellys04/projects/SmithLab-ChIPSeq_2016-12-31/analysis_dir/project_notes"




# ~~~~ CUSTOM FUNCITONs ~~~~~~ #
# make a custom track
def ucsc_cust_track_bigWig(track_name, track_url, view_settings):
	# track_template='track type=bigWig name="{}" bigDataUrl={} {}'.format(track_name, track_url, view_settings) # this doesn't work in python 2.6.6: 
	track_template='track type=bigWig name="%s" bigDataUrl=%s%s %s' % (track_name, track_url, track_name, view_settings)
	return track_template

# get the desired files from the directory
def get_files(dir_path,pattern):
	os.chdir(dir_path)
	got_files=[file for file in glob.glob(pattern)]
	return got_files

# write tracks for all input files to the output file (TSV format)
# def make_all_tracks(dir_path,pattern,track_url,view_settings,output_file):
	
	




# ~~~~ GET FILES ~~~~~~ #
# find the bigwig files in the external dir
chipseq_files = get_files(chipseq_external_bigwig_dir,"*.bw")

rnaseq_files = get_files(rnaseq_external_bigwig_dir,"*.bw")





# ~~~~ MAKE TRACKS ~~~~~~ #
# for long table single column of tracks
track_table_long = pd.DataFrame()
track_table_long["col1"] = pd.Series([ucsc_cust_track_bigWig(file,chipseq_external_URL,chipseq_bigwig_view_settings) for file in chipseq_files])
# write to the file
track_table_long.sort(["col1"]).to_csv(analysis_outdir + '/UCSC_ChIPSeq_tracks-long.tsv',sep='\t',index=False,header=False, quoting=csv.QUOTE_NONE)

track_table_long = pd.DataFrame()
track_table_long["col1"] = pd.Series([ucsc_cust_track_bigWig(file,rnaseq_external_URL,rnaseq_bigwig_view_settings) for file in rnaseq_files])
# write to the file
track_table_long.sort(["col1"]).to_csv(analysis_outdir + '/UCSC_RNASeq_tracks-long.tsv',sep='\t',index=False,header=False, quoting=csv.QUOTE_NONE)


quit()


# NOTES: 
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


# TO-DO: add more types of tracks! see the bash version of this for functionality to port over to Python
