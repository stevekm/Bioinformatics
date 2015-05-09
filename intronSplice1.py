from __future__ import division
DNA="ATCGATCGATCGATCGACTGACTAGTCATAGCTATGCATGTAGCTACTCGATCGATCGATCGATCGATCGATCGATCGATCGATCATGCTATCATCGATCGATATCGATGCATCGACTACTAT"
# The first exon runs from the start of the sequence to the sixty-third character, and the second exon runs from the ninety-first character to the end of the sequence. Write a program that will print just the coding regions of the DNA sequence.
# round off that percentage and display % symbol
DNALen=len(DNA)
print("DNA sequence is \n" + DNA)
print("Sequence length is " + str(DNALen))
print("Exon 1: " + DNA[0:62])
print("Exon 2: " + DNA[90:DNALen])
# write a program that will calculate what percentage of the DNA sequence is coding.
CodingRegionLength=len(DNA[0:62])+len(DNA[90:DNALen])
print("Coding region lenght is " + str(CodingRegionLength))
PercentCoding=CodingRegionLength/DNALen*100
print("Percent coding DNA is " + str(PercentCoding))
# write a program that will print out the original genomic DNA sequence with coding bases in uppercase and non-coding bases in lowercase.
print(DNA[0:62].upper() + DNA[63:89].lower() + DNA[90:DNALen].upper())
