#!/usr/bin/env python
# python 2.7

# "Sort a fasta file with ~5k sequences and match to a cluster txt file with sequence identifiers"
# output a fasta file that is sorted based on a list of ID's given in another file
# https://www.biostars.org/p/213003/
# http://biopython.org/wiki/SeqIO
# http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc16

from Bio import SeqIO

# pattern = sys.argv[1]
# file_path = sys.argv[2]
# fasta_file = "/local/data/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"


# $ head cluster.txt
# >Cluster 0
# >0  214aa, >gi|953781338|gb|ALO92858.1|... at 93.93%
# >1  211aa, >gi|504484658|ref|WP_014671760.1|... at 95.26%
# >2  211aa, >gi|664069898|ref|WP_030608495.1|... at 94.79%
# >3  213aa, >gi|927897990|ref|WP_053850251.1|... at 92.49%
# >4  213aa, >gi|973382757|ref|WP_059125750.1|... at 90.61%
# >5  211aa, >gi|664112575|ref|WP_030650029.1|... at 94.79%
# >6  210aa, >gi|518972364|ref|WP_020128239.1|... at 91.43%
# >7  217aa, >gi|822881081|emb|CQR64871.1|... at 91.71%

# # make a file with just the gi ID's as fasta format record ID's
# $ cat cluster.txt | sed -e 's/^>.*\(>.*\)$/\1/g' -e 's/^\(.*|\)\(\.\.\..*\)$/\1/'| grep -E -v '^>Cluster.*' > cluster_new.txt

# $ cat cluster_new.txt
# >gi|953781338|gb|ALO92858.1|
# >gi|504484658|ref|WP_014671760.1|
# >gi|664069898|ref|WP_030608495.1|
# >gi|927897990|ref|WP_053850251.1|
# >gi|973382757|ref|WP_059125750.1|
# >gi|664112575|ref|WP_030650029.1|
# >gi|518972364|ref|WP_020128239.1|
# >gi|822881081|emb|CQR64871.1|


# the paths to the input and output files
fasta_file = "/ifs/home/kellys04/biopython_tests/test.fa"
cluster_file = "/ifs/home/kellys04/biopython_tests/cluster_new.txt"
output_file = "/ifs/home/kellys04/biopython_tests/sorted.fa"


# start an empty list for the ID
ID_list = []

# make a list of the ID's
for record in SeqIO.parse(open(cluster_file, "rU"), "fasta"):
    # add each record ID to the list
    ID_list.append(record.id)

print ID_list


# iterate over the ID's, save the matching fasta entries
for my_ID in ID_list:
    # print my_ID
    for record in SeqIO.parse(open(fasta_file, "rU"), "fasta"):
        # if my_ID in record.id:
        if my_ID in record.id:
            # for appending to the output file
            output_handle = open(output_file, "a")
            SeqIO.write(record, output_handle, "fasta")
            output_handle.close()



