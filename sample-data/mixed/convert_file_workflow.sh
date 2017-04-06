#!/bin/bash

chrom_sizes="/ifs/home/kellys04/software/hg19.chrom.sizes"

# start with test.bam

# sort
samtools sort test.bam -o test.bam.sorted
/bin/mv test.bam.sorted test.bam

# index
samtools index test.bam

# convert to bed
samtools view test.bam | awk '{print $3 "\t" $4 "\t" $4+length($10)-1}' > test.bed
# bamToBed -i test.bam > test.bed

# convert to bigBed
bedToBigBed test.bed "$chrom_sizes" test.bigBed

# convert to bedGraph
bedtools genomecov -ibam test.bam -g "$chrom_sizes" -bg > test.bedgraph
# genomeCoverageBed –i file.bed -bg –g my.genome > sample.cov

# convert to bigWig
bedGraphToBigWig test.bedgraph "$chrom_sizes" test.bigwig
