#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Subset a series of fastq.gz files based on reads that fall in desired target regions
use a .bam file of reads from the fastq's and a .bed file of target regions
to find which reads fall within the desired regions
then print those same reads from the original fastq to a new fastq file


http://biopython.org/wiki/SeqIO
http://pysam.readthedocs.io/en/latest/api.html

"""
import os
import sys
import time
import csv
import gzip
from Bio import SeqIO
import pysam

# ~~~~~ CLASSES ~~~~~ #
class TimeTracker(object):
    """
    Object to track time in the program

    Examples
    --------
    Example usage::

        time_tracker = TimeTracker()
        time_tracker.elapsed()
    """
    def __init__(self, print_start = False):
        self.start_time = time.time()
        if print_start:
            self.started()
    def started(self):
        """
        Prints the time when the object was started
        """
        print('start time: {0}'.format(self.start_time)); sys.stdout.flush() # flush messages immediately or they hang in the buffer a long time
    def elapsed(self):
        """
        Sets and prints the elapsed time since the obect creation
        """
        self.started()
        self.end_time = time.time()
        print('current time: {0}'.format(self.end_time)); sys.stdout.flush()
        self.time_elapsed = self.end_time - self.start_time
        print('time_elapsed: {0} seconds ({1} minutes)'.format(self.time_elapsed, int(self.time_elapsed) / 60)); sys.stdout.flush()


# ~~~~~ FUNCTIONS ~~~~~ #
def read_targets(target_file):
    """
    Reads the genomic regions from a .bed file

    Parameters
    ----------
    target_file: str
        path to a tab-separated .bed formatted file of genomic coordinates

    Returns
    -------
    list
        a list of tuples of the format 'chrom', 'start', 'stop'
    """
    # read in the target regions
    targets = []
    with open(target_file) as f:
        for line in f:
            chrom, start, stop = line.split()
            targets.append((chrom, int(start), int(stop)))
    return(targets)

def bam_target_qnames(bam_file, targets):
    """
    Retrieves the `qname` sequence ID's for reads in a .bam file that are at specified target coordinates

    Parameters
    ----------
    bam_file: str
        path to a .bam file
    targets:
        a list of tuples with 'chrom', 'start', and 'stop' coordinates

    Returns
    -------
    list
        a list of qnames IDs
    """
    qnames = []
    for chrom, start, stop in targets:
        # load the bam
        bam = pysam.AlignmentFile(bam_file, "rb")
        for read in bam.fetch(chrom, start, stop):
            qnames.append(read.qname)
        bam.close()
    return(qnames)

def get_all_bam_qnames(bam_files, targets):
    """
    Gets all the qnames for a set of bam files

    Parameters
    ----------
    bam_files: list
        list of pathes to .bam files
    targets:
        a list of tuples with 'chrom', 'start', and 'stop' coordinates

    Returns
    -------
    list
        a list of qnames IDs
    """
    qnames = []
    for bam_file in bam_files:
        print('reading from .bam file: {0}'.format(bam_file))
        for q in bam_target_qnames(bam_file = bam_file, targets = targets):
            qnames.append(q)
    qnames = sorted(set(qnames))
    print('read {0} unique qnames'.format(len(qnames)))
    return(qnames)


def filter_fastq_qnames(fastq_in, fastq_out, qnames):
    """
    Outputs entries in a .fastq.gz file that match the qname sequence ID into a new fastq.gz file

    Parameters
    ----------
    fastq_in: str
        path to .fastq.gz file to read sequences from
    fastq_out: str
        path to .fastq.gz file to write sequences to
    qnames: list
        list of qname ID's to search for in the .fastq
    """
    print('Fastq in: {0}\nFastq out: {1}\n\n'.format(fastq_in, fastq_out)); sys.stdout.flush()
    with gzip.open(fastq_in) as gz_in, gzip.open(fastq_out, 'wb') as gz_out:
        input_seq_iterator = SeqIO.parse(gz_in, "fastq")
        seq_iterator = (record for record in input_seq_iterator if record.id in qnames)
        SeqIO.write(seq_iterator, gz_out, "fastq")

def filter_all_fastqs(fastqs, qnames, output_dir = 'output_fastq'):
    """
    Filters all the fastq files for entries in a that match the qname sequence IDs

    Parameters
    ----------
    fastqs: list
        a list of paths to .fastq.gz files to be filtered
    output_dir: str
        path to output directory
    qnames: list
        list of qname ID's to search for in the .fastq
    """
    for fastq in fastqs:
        fastq_out = os.path.join(output_dir, os.path.basename(fastq))
        filter_fastq_qnames(fastq_in = fastq, fastq_out = fastq_out, qnames = qnames)


def main():
    """
    Main control funciton for the script

    Examples
    --------
    Example usage::

        ./run_script.py input_fastq/HapMap-B17-1267_S8_L002_R2_001.fastq.gz,input_fastq/HapMap-B17-1267_S8_L003_R1_001.fastq.gz,input_fastq/HapMap-B17-1267_S8_L002_R1_001.fastq.gz,input_fastq/HapMap-B17-1267_S8_L001_R1_001.fastq.gz,input_fastq/HapMap-B17-1267_S8_L001_R2_001.fastq.gz,input_fastq/HapMap-B17-1267_S8_L004_R2_001.fastq.gz,input_fastq/HapMap-B17-1267_S8_L004_R1_001.fastq.gz,input_fastq/HapMap-B17-1267_S8_L003_R2_001.fastq.gz input_bam/HapMap-B17-1267.bam

    """
    # file with genomic regions to use
    target_file = "targets.bed"
    # directory to save the output fastq files
    output_dir= "output_fastq"


    args = sys.argv[1:]

    fastqs = args[0].split(',')

    bams = args[1].split(',')


    sample = {
    'fastq': fastqs, # list of fastq file paths
    'bam': bams, # list of bam file paths
    'qname': [] # list to hold seq ID's from the bam to look up in the fastq
    }

    print(sample); sys.stdout.flush()

    # get the targets
    targets = read_targets(target_file = target_file)

    # get the qnames
    sample['qname'] = get_all_bam_qnames(bam_files = sample['bam'], targets = targets)

    # filter the fastq files
    filter_all_fastqs(fastqs = sample['fastq'], qnames = sample['qname'], output_dir = output_dir)


# ~~~~~ RUN ~~~~~ #
if __name__ == '__main__':
    time_tracker = TimeTracker()
    main()
    time_tracker.elapsed()
