# this program generates the complimentary sequence of the given DNA sequence
# in future versions allow user input of DNA sequence
DNA="ACTGATCGATTACGTATAGTATTTGCTATCATACATATATATCGATGCGTTCAT"
print(DNA)
print("complimentary sequence")
cDNA=DNA.replace("A","%temp%").replace("T","A").replace("%temp%","T").replace("G","%temp%").replace("C","G").replace("%temp%","C")
print(cDNA)