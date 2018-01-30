# makefile to download iGenomes data for hg19 
# https://support.illumina.com/sequencing/sequencing_software/igenome.html
SHELL:=/bin/bash
GENOME_FA_MD5:=c1ddcc5db31b657d167bea6d9ff354f9

none:

Homo_sapiens_UCSC_hg19.tar.gz:
	wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg19/Homo_sapiens_UCSC_hg19.tar.gz

list: 
	wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg19/Homo_sapiens_UCSC_hg19.tar.gz -O- | tar -ztvf - > tar_contents.txt

genome.fa:
	wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg19/Homo_sapiens_UCSC_hg19.tar.gz -O- | tar -ztvf - Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
