#!/bin/bash

# PROBLEM: want to preseve the terminal output from bowtie2, but also copy the stderr from bowtie into a separate file
# SOLUTION: use `tee` along with some bash stream redirection to copy the stderr stream to a new file AND print it on the terminal


# set files and places
tmp_fastq1="$HOME/projects/SmithLab_PARCLIP/14Q-sample1_R1_mini.fastq"
tmp_fastq2="$HOME/projects/SmithLab_PARCLIP/14Q-sample1_R2_mini.fastq"
tmp_outdir="$HOME/projects/SmithLab_PARCLIP/test_bowtie2"
cd $tmp_outdir
module load bowtie2 # http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
module load samtools/1.2.1
mapq=30
THREADS=${NSLOTS:=8}
tmpGenome="/local/data/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"
tmpOUTFILE="$(basename "$tmp_fastq1")"

# ~~~~~~~~~ # 
# FINAL COMMAND TO USE:
(bowtie2 --threads "$THREADS" --local -x "$tmpGenome" -q -1 "$tmp_fastq1" -2 "$tmp_fastq2" --met-file bowtie2_alignment-metrics.txt | samtools view -@ "$THREADS" -Sb1 - | samtools sort -m 10G -@ "$THREADS" - "$tmpOUTFILE" ) 3>&1 1>&2 2>&3 | tee stderr.log 

# ~~~~~~~~~ #
# EXPLANATIONS & EXAMPLES:
# standard bowtie2 alignment command, pipes output directly to samtools for sorting and conversion to .bam file; stderr is printed on screen and contains alignment stats from bowtie; no stderr or stdout from samtools
bowtie2 --threads "$THREADS" --local -x "$tmpGenome" -q -1 "$tmp_fastq1" -2 "$tmp_fastq2" --met-file bowtie2_alignment-stats.txt | samtools view -@ "$THREADS" -Sb1 - | samtools sort -m 10G -@ "$THREADS" - "$tmpOUTFILE"


# redirect just the stderr from bowtie to a file
(bowtie2 --threads "$THREADS" --local -x "$tmpGenome" -q -1 "$tmp_fastq1" -2 "$tmp_fastq2" --met-file bowtie2_alignment-metrics.txt) 2>stats.txt


# redirect bowtie stderr to text file while piping stdout
# this does not give console output though, which we want to keep
(bowtie2 --threads "$THREADS" --local -x "$tmpGenome" -q -1 "$tmp_fastq1" -2 "$tmp_fastq2" --met-file bowtie2_alignment-metrics.txt) 2>stats.txt | samtools view -@ "$THREADS" -Sb1 - | samtools sort -m 10G -@ "$THREADS" - "$tmpOUTFILE"


# without stderr redirection; stderr prints to console
bowtie2 --threads "$THREADS" --local -x "$tmpGenome" -q -1 "$tmp_fastq1" -2 "$tmp_fastq2" --met-file bowtie2_alignment-metrics.txt | samtools view -@ "$THREADS" -Sb1 - | samtools sort -m 10G -@ "$THREADS" - "$tmpOUTFILE"


# copy stderr while piping stdout
# (cmd | tee stdout.log) 3>&1 1>&2 2>&3 | tee stderr.log
# # copy stdout to a log file 
# # (tee only accepts stdout as input!)
(bowtie2 --threads "$THREADS" --local -x "$tmpGenome" -q -1 "$tmp_fastq1" -2 "$tmp_fastq2" --met-file bowtie2_alignment-metrics.txt | tee stdout.log )


# # copy stdout AND stderr to stderr to separate log files
# 3>&1 1>&2 2>&3
# # redirect '3' to stdout, redirect stdout ('1') to stder ('2'), redirect stderr to '3'
(bowtie2 --threads "$THREADS" --local -x "$tmpGenome" -q -1 "$tmp_fastq1" -2 "$tmp_fastq2" --met-file bowtie2_alignment-metrics.txt | tee stdout.log ) 3>&1 1>&2 2>&3 | tee stderr.log 
# # first bowtie stdout is getting copied to stdout.log with tee 
# # then the new stream '3' is getting redirected to stdout, stderr gets redirected to '3' which nows goes to stdout (for tee), and stdout gets sent to stderr
# # the stdout stream, which now contains only the original stderr stream, gets copied into stderr.log with tee; the new stdout stream now only contains original stderr
# # '3' gets lost after the pipe


# you can see that the final tee's incoming stdout stream is the original stderr stream
(bowtie2 --threads "$THREADS" --local -x "$tmpGenome" -q -1 "$tmp_fastq1" -2 "$tmp_fastq2" --met-file bowtie2_alignment-metrics.txt | tee stdout.log ) 3>&1 1>&2 2>&3 | tee stderr.log > tee_tmp.txt


# the '3' stream no longer exists after the pipe
(bowtie2 --threads "$THREADS" --local -x "$tmpGenome" -q -1 "$tmp_fastq1" -2 "$tmp_fastq2" --met-file bowtie2_alignment-metrics.txt | tee stdout.log ) 3>&1 1>&2 2>&3 | tee stderr.log 3> tee_tmp.txt


# THIS IS THE ONE WE WANT! FINAL COMMAND TO USE:
# redirect just the stderr to a log file, don't need to keep stdout from bowtie
(bowtie2 --threads "$THREADS" --local -x "$tmpGenome" -q -1 "$tmp_fastq1" -2 "$tmp_fastq2" --met-file bowtie2_alignment-metrics.txt | samtools view -@ "$THREADS" -Sb1 - | samtools sort -m 10G -@ "$THREADS" - "$tmpOUTFILE" ) 3>&1 1>&2 2>&3 | tee stderr.log 


# if you want to redirect both the stderr and the stdout to logs (warning: stdout.log will contain the SAM file, huge!)
(bowtie2 --threads "$THREADS" --local -x "$tmpGenome" -q -1 "$tmp_fastq1" -2 "$tmp_fastq2" --met-file bowtie2_alignment-metrics.txt | tee stdout.log | samtools view -@ "$THREADS" -Sb1 - | samtools sort -m 10G -@ "$THREADS" - "$tmpOUTFILE" ) 3>&1 1>&2 2>&3 | tee stderr.log 


# RESOURCES:
# http://stackoverflow.com/questions/692000/how-do-i-write-stderr-to-a-file-while-using-tee-with-a-pipe
# http://unix.stackexchange.com/questions/6430/how-to-redirect-stderr-out-to-different-files-and-also-display-in-terminal

