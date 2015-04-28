from __future__ import division
# this program repoorts the nucleotide content of a given DNA sequence
# in future versions, allow for user input of DNA sequence
# in future versions, allow for user input of desired nucleotide(s) for analysis
# in future versions, round off output and display '%' symbol
DNA="ACTGATCGATTACGTATAGTATTTGCTATCATACATATATATCGATGCGTTCAT"
DNAlength=len(DNA)
#print(len(DNA))
#print(DNAlength)
print("Sequence Length = " + str(DNAlength))
Acontent=DNA.count("A")
Tcontent=DNA.count("T")
ATcontent=(Acontent+Tcontent)/DNAlength*100
print("The number of A is " + str(Acontent))
print("The number of T is " + str(Tcontent))
print("The AT Content is " + str(ATcontent))
