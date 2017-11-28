#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Subset a series of fastq.gz files based on reads that fall in desired target regions
use a .bam file of reads from the fastq's and a .bed file of target regions
to find which reads fall within the desired regions
then print those same reads from the original fastq to a new fastq file


http://biopython.org/wiki/SeqIO
http://pysam.readthedocs.io/en/latest/api.html

$ cat sample_files.tsv
sample	bam	fastq
HapMap-B17-1267	input_bam/HapMap-B17-1267.bam	input_fastq/HapMap-B17-1267_S8_L002_R2_001.fastq.gz,input_fastq/HapMap-B17-1267_S8_L003_R1_001.fastq.gz,input_fastq/HapMap-B17-1267_S8_L002_R1_001.fastq.gz,input_fastq/HapMap-B17-1267_S8_L001_R1_001.fastq.gz,input_fastq/HapMap-B17-1267_S8_L001_R2_001.fastq.gz,input_fastq/HapMap-B17-1267_S8_L004_R2_001.fastq.gz,input_fastq/HapMap-B17-1267_S8_L004_R1_001.fastq.gz,input_fastq/HapMap-B17-1267_S8_L003_R2_001.fastq.gz
NTC-H2O	input_bam/NTC-H2O.bam	input_fastq/NTC-H2O_S1_L001_R2_001.fastq.gz,input_fastq/NTC-H2O_S1_L004_R2_001.fastq.gz,input_fastq/NTC-H2O_S1_L002_R1_001.fastq.gz,input_fastq/NTC-H2O_S1_L003_R2_001.fastq.gz,input_fastq/NTC-H2O_S1_L003_R1_001.fastq.gz,input_fastq/NTC-H2O_S1_L004_R1_001.fastq.gz,input_fastq/NTC-H2O_S1_L002_R2_001.fastq.gz,input_fastq/NTC-H2O_S1_L001_R1_001.fastq.gz
SeraCare-1to1-Positive	input_bam/SeraCare-1to1-Positive.bam	input_fastq/SeraCare-1to1-Positive_S2_L004_R1_001.fastq.gz,input_fastq/SeraCare-1to1-Positive_S2_L003_R2_001.fastq.gz,input_fastq/SeraCare-1to1-Positive_S2_L001_R1_001.fastq.gz,input_fastq/SeraCare-1to1-Positive_S2_L004_R2_001.fastq.gz,input_fastq/SeraCare-1to1-Positive_S2_L002_R2_001.fastq.gz,input_fastq/SeraCare-1to1-Positive_S2_L002_R1_001.fastq.gz,input_fastq/SeraCare-1to1-Positive_S2_L001_R2_001.fastq.gz,input_fastq/SeraCare-1to1-Positive_S2_L003_R1_001.fastq.gz
"""
import time
import os
import sys
import csv
import gzip
from Bio import SeqIO
import pysam


start_time = time.time()
print('start time: {0}'.format(start_time)); sys.stdout.flush() # flush messages immediately or they hang in the buffer a long time

# file with sample IDs, paths to .fastq and .bam files
samples_file = "sample_files.tsv"
# file with genomic regions to use
target_file = "targets.bed"
# directory to save the output fastq files
output_dir= "output_fastq"


# read in the samples
samples = []
with open(samples_file) as f:
    reader = csv.DictReader(f, delimiter = '\t')
    for sample in reader:
        samples.append(sample)

# split the fastq and bam entries
for sample in samples:
    sample['fastq'] = sample['fastq'].split(',')
    sample['bam'] = sample['bam'].split(',')
    # add 'qname' for later
    sample['qname'] = []
    # make fastq output file paths
    sample['output_fastq'] = {}
    for fastq in sample['fastq']:
        sample['output_fastq'][fastq] = os.path.join(output_dir, os.path.basename(fastq)) # .strip('.gz')

print(samples); sys.stdout.flush()



# read in the target regions
targets = []
with open(target_file) as f:
    for line in f:
        chrom, start, stop = line.split()
        targets.append((chrom, int(start), int(stop)))


# get the read IDs ('qname') for the targets from each bam per sample
for chrom, start, stop in targets:
    print((chrom, start, stop)); sys.stdout.flush()
    for sample in samples:
        for bam_file in sample['bam']:
            print(bam_file); sys.stdout.flush()
            # load the bam
            bam = pysam.AlignmentFile(bam_file, "rb")
            for read in bam.fetch(chrom, start, stop):
                sample['qname'].append(read.qname)
            bam.close()


# get only the unique qnames per sample
for sample in samples:
    sample['qname'] = sorted(set(sample['qname']))

# get the fastq reads that match
for sample in samples:
    for fastq in sample['fastq']:
        output_fastq = sample['output_fastq'][fastq]
        print((fastq, output_fastq)); sys.stdout.flush()
        with gzip.open(fastq) as gz_in, gzip.open(output_fastq, 'wb') as gz_out:
            input_seq_iterator = SeqIO.parse(gz_in, "fastq")
            seq_iterator = (record for record in input_seq_iterator if record.id in sample['qname'])
            SeqIO.write(seq_iterator, gz_out, "fastq")


end_time = time.time()
print('end_time: {0}'.format(end_time)); sys.stdout.flush()

time_elapsed = end_time - start_time
print('time_elapsed: {0}'.format(time_elapsed)); sys.stdout.flush()
