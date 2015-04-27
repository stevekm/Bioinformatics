from __future__ import division
# this program calculates the AT content of a pre-defined DNA sequence
# in future versions, allow for user input of DNA sequence
# round off the AT content and display '%'
DNA="ACTGATCGATTACGTATAGTATTTGCTATCATACATATATATCGATGCGTTCAT"
print("DNA sequence is " + DNA)
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
