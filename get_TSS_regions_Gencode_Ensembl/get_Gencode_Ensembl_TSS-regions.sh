#!/bin/bash
set -x 
# this is the pipeline used to get the TSS regions from the Gencode and Ensembl reference files


# ~~~~~ SOFTWARE ~~~~~~~# 
# load some programs to be used; this updates the PATH with the corresponding entries
module load homer/v4.6
module unload gcc
module load bedtools/2.22.0

# install BEDOPS http://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/gtf2bed.html
# wget https://github.com/bedops/bedops/releases/download/v2.4.16/bedops_linux_x86_64-v2.4.16.tar.bz2; tar jxvf bedops_linux_x86_64-vx.y.z.tar.bz2
# all the BEDOPS binaries are unpacked in this location, update the PATH for it:
PATH=$PATH:~/software/bin

# use gtools to get the TSS regions
# I think this is it? https://github.com/tsirigos/ibm-cbc-genomic-tools/blob/master/gtools/genomic_regions.cpp
# VERSION: genomic-tools 3.0.0
# this is already in my PATH

# ~~~~~ COMMON ITEMS ~~~~~~~# 
# dir for outputs and such
ProjDir="~/Gencode_Ensembl_TSS_regions"
mkdir -p "$ProjDir"; cd "$ProjDir"

# size of the TSS regions to use
region_size="10000"



# ~~~~~ GENCODE STEPS ~~~~~~~# 
# file prefix
gen_file="gencode.v19.annotation"

# Download source data ; http://www.gencodegenes.org/releases/19.html
[ ! -f ${gen_file}.gtf.gz ] && wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz

# convert the GTF to BED witht he BEDOPS converter
[ ! -f ${gen_file}.bed ] && zcat ${gen_file}.gtf.gz | convert2bed --input=gtf - > ${gen_file}.bed

# use gtools to get the TSS regions
# VERSION: genomic-tools 3.0.0
[ ! -f ${gen_file}.tss.bed ] && cat ${gen_file}.bed | genomic_regions connect | genomic_regions pos -op 5p > ${gen_file}.tss.bed

# collapse the peaks with HOMER
[ ! -f ${gen_file}.tss_collapse.txt ] && mergePeaks ${gen_file}.tss.bed -strand > ${gen_file}.tss_collapse.txt
[ ! -f ${gen_file}.tss_collapse.bed ] && pos2bed.pl ${gen_file}.tss_collapse.txt > ${gen_file}.tss_collapse.bed && rm -f ${gen_file}.tss_collapse.txt

# create 10,000kbp TSS regions with bedtools
if [ ! -f ${gen_file}.tss_collapse_${region_size}.bed ]; then
  cat ${gen_file}.tss_collapse.bed | bedtools slop -g hg19.chrom.sizes -b "$region_size" > ${gen_file}.tss_collapse_${region_size}.bed
fi

# ~~~~~~~~~~ # ~~~~~~~~~~ # ~~~~~~~~~~ # ~~~~~~~~~~ # ~~~~~~~~~~ 


# ~~~~~ ENSEMBL STEPS ~~~~~ #
# file prefix
ens_file="genes.ensembl.GRCh37.82"
# download source data
if [ ! -f ${ens_file}.gtf.gz ]; then
  wget ftp://ftp.ensembl.org/pub/grch37/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh37.82.gtf.gz -O ${ens_file}.gtf.gz
fi

# GTF processing
# Need to filter out the GL, MT entries, need to add 'chr' to the start of each 1st entry
if [ ! -f ${ens_file}_noGLMT.gtf ]; then
  zcat ${ens_file}.gtf.gz | grep -Ev "^#|^GL|^M" > ${ens_file}_noGLMT.gtf
  # cat ${ens_file}_noGLMT.gtf | awk '{ OFS="\t"; $1 = "chr"$1; print }' > ${ens_file}_noGLMT_chr.gtf # this don't work?? 
  cat ${ens_file}_noGLMT.gtf | sed 's/^/chr/' > ${ens_file}_noGLMT_chr.gtf
  
fi

# convert to BED 
if [ ! -f ${ens_file}_noGLMT_chr.bed ]; then
  # cat ${ens_file}_noGLMT_chr.gtf | convert2bed --input=gtf - > ${ens_file}_noGLMT_chr.bed
  gtf2bed < ${ens_file}_noGLMT_chr.gtf > ${ens_file}_noGLMT_chr.bed
fi

# get the TSS regions
if [ ! -f ${ens_file}_noGLMT_chr.tss.bed ]; then
  cat ${ens_file}_noGLMT_chr.bed | genomic_regions connect | genomic_regions pos -op 5p > ${ens_file}_noGLMT_chr.tss.bed
fi

# collapse the regions
 if [ ! -f ${ens_file}_noGLMT_chr.tss_collapse.bed ]; then
  mergePeaks ${ens_file}_noGLMT_chr.tss.bed -strand > ${ens_file}_noGLMT_chr.tss_collapse.txt
  pos2bed.pl ${ens_file}_noGLMT_chr.tss_collapse.txt > ${ens_file}_noGLMT_chr.tss_collapse.bed && rm -f ${ens_file}_noGLMT_chr.tss_collapse.txt
fi

# create 10,000kbp TSS regions
if [ ! -f ${ens_file}_noGLMT_chr.tss_collapse_${region_size}.bed ]; then
  cat ${ens_file}_noGLMT_chr.tss_collapse.bed | bedtools slop -g hg19.chrom.sizes -b "$region_size" > ${ens_file}_noGLMT_chr.tss_collapse_${region_size}.bed
fi



# Here are the chromosome sizes for the reference hg19 genomes used

# $ cat hg19.chrom.sizes
# chrM	16571
# chr1	249250621
# chr2	243199373
# chr3	198022430
# chr4	191154276
# chr5	180915260
# chr6	171115067
# chr7	159138663
# chr8	146364022
# chr9	141213431
# chr10	135534747
# chr11	135006516
# chr12	133851895
# chr13	115169878
# chr14	107349540
# chr15	102531392
# chr16	90354753
# chr17	81195210
# chr18	78077248
# chr19	59128983
# chr20	63025520
# chr21	48129895
# chr22	51304566
# chrX	155270560
# chrY	59373566
