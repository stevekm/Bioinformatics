#!/usr/bin/python
from pybedtools.contrib.bigwig import bam_to_bigwig
import sys

# print 'Number of arguments:', len(sys.argv), 'arguments.', '\n'
# print 'Argument List:', str(sys.argv), '\n'
print 'Python script follows here...'
###

print 'Genome is: ', str(sys.argv[1]), '\n'

print 'Bam input file is: ', str(sys.argv[2]), '\n'

print 'Bigwig output file is: ', str(sys.argv[3]), '\n'


# bam_to_bigwig(bam='path/to/bam', genome='mm9', output='path/to/bigwig')
bam_to_bigwig(bam=str(sys.argv[2]), genome=str(sys.argv[1]), output=str(sys.argv[3]))

