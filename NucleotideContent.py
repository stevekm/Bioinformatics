from __future__ import division
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
