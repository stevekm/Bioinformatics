#!/bin/bash
set -x
##
## USAGE: HITSCLIP_pipeline_trim16_qsub.sh /path/to/fastqc_outdir /path/to/input_file /path/to/novoindex/genome <sampleID> /path/to/ref_genome.fa /path/to/meme_db /path/to/meme_db2
## This script will run the comlete HITS-CLIP pipeline on a single sample
## This script should be submitted with qsub to the SGE cluster, see the accompanying file
## 'workflow' for setup scripts to parse out the correct arguments per sample and pass them to this script, and setup
## the needed directory structure

# ~~~~~~ permissions ~~~~~~ #
# To make stuff group-writeable (this is what I add to all my shared scripts):
umask 007


# ~~~~~~ script args processing ~~~~~~ #
# # if not enough args, output USAGE info lines (which start with '##') and exit
if (($# != 7)); then
  grep '^##' $0
  exit
fi


# make sure that the correct number of script arguments were used
# if [ $# != 2 ] # if not enough args provided
# then
#   grep '^##' $0 # print out lines from the script that start with '##' and exit
#   exit
# fi


# get the arguments
OUTDIR="$1" # outdir
INPUTFILE="$2" # input file
GENOME="$3" # file path to ref genome for Novoalign
SAMPLEID="$4"
GENOMEFA="$5" # file path to ref genome fasta file
tmpMEME_DB="$6" # file path to the MEME DB to use
tmpMEME_DB2="$7"

THREADS=${NSLOTS:=8}
echo "Working dir is $PWD"
echo "OUTDIR is $OUTDIR"
echo "INPUTFILE is $INPUTFILE"
echo "GENOME is $GENOME"
echo "SAMPLEID is $SAMPLEID"
echo "GENOMEFA is $GENOMEFA"
echo "tmpMEME_DB is $tmpMEME_DB"
echo "tmpMEME_DB2 is $tmpMEME_DB2"
echo "THREADS is $THREADS"

# get the file name
tmp_filename=$(basename "$INPUTFILE" )
tmp_file_base=${tmp_filename%%.*}
echo "tmp_filename is $tmp_filename"
echo "tmp_file_base is $tmp_file_base"
# strip the extension
# tmp_samplename="${tmp_filename%%.fasta}" 


# software locations dirs
CIMSdir="/ifs/home/kellys04/software/CIMS"
NOVOdir="/ifs/home/kellys04/software/novocraft"

# load HOMER
module load homer/v4.6
module load fastqc/0.11.4
export KRAKEN_DEFAULT_DB="$HOME/ref/Kraken/nt-48G"
export KRAKEN_NUM_THREADS="$THREADS"


# ~~~~~~ Script Logging ~~~~~~~ #
# get the script file 
zz=$(basename $0)
# get the script dir
za=$(dirname $0)
# use either qsub job name or scriptname in log filename
mkdir -p logs/scripts
LOG_FILE=logs/scripts/scriptlog.${JOB_NAME:=$(basename $0)}.$(date -u +%Y%m%dt%H%M%S).${HOSTNAME:-$(hostname)}.$$.${RANDOM}
echo "This is the log file for the script." > $LOG_FILE
echo -e "\n pwd is:\n$PWD\n" >> $LOG_FILE
# echo -e "\nScript file is:\n$0\n" >> $LOG_FILE # for regular script usage (no qsub)
echo -e "\nScript file is:\n${JOB_NAME:=$(readlink -m $0)}\n" >> $LOG_FILE 
echo -e "\nScript file contents:\n\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%" >> $LOG_FILE
cat "$0" >> $LOG_FILE
# ~~~~~~~~~~~~ #

# ~~~~~~ # ~~~~~~ # ~~~~~~ run commands ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ 

# ~~~~~~ # ~~~~~~  QUALITY CONTROL # ~~~~~~ # ~~~~~~ # ~~~~~~ 
# FastQC
# FASTQC_DIR="${tmp_filename}_fastqc"
# if [ ! -d ${FASTQC_DIR} ]; then
#   mkdir -p "$FASTQC_DIR"
#   fastqc --threads "$THREADS" --nogroup --outdir "$FASTQC_DIR" "$INPUTFILE"
# fi
# # set a QC outdir
tmpQCdir="${OUTDIR}/quality_control"
if [ ! -d $tmpQCdir ]; then
  mkdir -p "$tmpQCdir"

  # FastQC; check raw data quality
  tmpFastQCdir1="${tmpQCdir}/FastQC_$(basename $INPUTFILE)"
  mkdir -p "$tmpFastQCdir1"
  fastqc --threads "$THREADS" --nogroup --outdir "$tmpFastQCdir1" "$INPUTFILE"
fi



# ~~~~~~ # ~~~~~~ HITS-CLIP Custom Pipeline # ~~~~~~ # ~~~~~~ # ~~~~~~ # ~~~~~~ 
# filter out low quality reads
if [ ! -f ${SAMPLEID}.fa ]; then
  echo -e "Filtering low quality reads"
  perl $CIMSdir/fastq_filter.pl -v -f mean:0-24:20 -of fasta "$INPUTFILE" "${SAMPLEID}.fa"
fi

# collapse duplicates
if [ ! -f ${SAMPLEID}.c.fa ]; then
  echo -e "Collapsing duplicate reads"
  perl $CIMSdir/fasta2collapse.pl -v "${SAMPLEID}.fa" "${SAMPLEID}.c.fa"
fi

# # trim back the first 20 nt because they are bad
if [ ! -f ${SAMPLEID}.c.t16.fa ]; then
  echo -e "Trimming reads"
  /ifs/home/kellys04/software/bin/fastx_trimmer -f 16 -i "${SAMPLEID}.c.fa" -o "${SAMPLEID}.c.t16.fa"
fi
 

# Kraken; check for potential contaminants
tmpKrakendir="${tmpQCdir}/Kraken"
# # use this for running on raw fastq.gz files
# tmp_R1="${tmpKrakendir}/$(basename $INPUTFILE)"
# # # only use the first 1,000,000 reads !!
# zcat $INPUTFILE | head -4000000 | gzip > "$tmp_R1"
# $HOME/software/bin/kraken --gzip-compressed --fastq-input $tmp_R1 | $HOME/software/bin/kraken-report | awk -F $'\t' '$1>0.1' > "${tmpKrakendir}/kraken_paired.$(basename $INPUTFILE).txt"
if [ ! -f ${tmpKrakendir}/${SAMPLEID}_kraken.txt ]; then
  mkdir -p "$tmpKrakendir"
  $HOME/software/bin/kraken --fasta-input  "${SAMPLEID}.c.t16.fa" | $HOME/software/bin/kraken-report | awk -F $'\t' '$1>0.1' > "${tmpKrakendir}/${SAMPLEID}_kraken.txt"
fi



# align with NovoAlign
# check if the file exists already because this takes a while
if [ ! -f ${SAMPLEID}.novoalign ]; then
  echo -e "Aligning reads"
  $NOVOdir/novoalign -t 85 -d "$GENOME" -f  "${SAMPLEID}.c.t16.fa" -F FA -l 25 -s 1 -r None -c "$THREADS" > "${SAMPLEID}.novoalign"
fi


# 
# Parse the ‘novoalign’ output and save coordinates of unambiguously mapped tags and mutations detected in these tags 
# in two separate files:
if [ ! -f ${SAMPLEID}.novoalign.tag.bed ]; then
  echo -e "Saving coordinates of single mapped alignments, mutations"
  perl $CIMSdir/novoalign2bed.pl -v --mismatch-file "${SAMPLEID}.novoalign.mutation.txt" "${SAMPLEID}.novoalign" "${SAMPLEID}.novoalign.tag.bed"
fi

# Collapse the potential PCR duplicates by coordinates and identify unique CLIP tags. 
if [ ! -f ${SAMPLEID}.novoalign.tag.uniq.bed ]; then
  echo -e "Collapsing PCR duplicates"
  perl $CIMSdir/tag2collapse.pl -v -weight --weight-in-name --keep-max-score --keep-tag-name "${SAMPLEID}.novoalign.tag.bed" "${SAMPLEID}.novoalign.tag.uniq.bed"
fi

#  Prepare a bedGraph file of unique tags for visualization in genome browser
if [ ! -f ${SAMPLEID}.novoalign.tag.uniq.bedgraph ]; then
  echo -e "Creating bedGraph"
  start_time1=`date +%s`
  perl $CIMSdir/tag2profile.pl -v -ss -exact -of bedgraph -n "Unique Tag Profile" "${SAMPLEID}.novoalign.tag.uniq.bed" "${SAMPLEID}.novoalign.tag.uniq.bedgraph"
  echo -e "Execution time:\t$(expr $(date +%s) - $start_time1)s"
fi

#  Cluster overlapping CLIP tags:
if [ ! -f ${SAMPLEID}.novoalign.tag.uniq.cluster.bed ]; then
  echo -e "Clustering the overlapping CLIP tags"
  start_time1=`date +%s`
  perl $CIMSdir/tag2cluster.pl -v -s -maxgap "-1" "${SAMPLEID}.novoalign.tag.uniq.bed" "${SAMPLEID}.novoalign.tag.uniq.cluster.bed"
  echo -e "Execution time:\t$(expr $(date +%s) - $start_time1)s"
fi

# 
#  Extract PH of each cluster for ranking:
# The fifth column of the output file represents the PH in each cluster.
if [ ! -f ${SAMPLEID}.novoalign.tag.uniq.cluster.PH.bed ]; then
  echo -e "Extracting PH of each cluster" # what is a PH ????
  perl $CIMSdir/extractPeak.pl -s -v "${SAMPLEID}.novoalign.tag.uniq.cluster.bed" "${SAMPLEID}.novoalign.tag.uniq.bedgraph" "${SAMPLEID}.novoalign.tag.uniq.cluster.PH.bed"
fi

#  Extract mutations in unique CLIP tags 
if [ ! -f ${SAMPLEID}.novoalign.tag.uniq.mutation.txt ]; then
  echo -e "Extracting unique mutations"
  python $CIMSdir/joinWrapper.py ${SAMPLEID}.novoalign.mutation.txt ${SAMPLEID}.novoalign.tag.uniq.bed 4 4 N ${SAMPLEID}.novoalign.tag.uniq.mutation.txt
fi

echo "Separating different types of mutations:"
if [ ! -f ${SAMPLEID}.novoalign.tag.uniq.del.bed ]; then awk '{if($9=="-") {print $0}}' ${SAMPLEID}.novoalign.tag.uniq.mutation.txt | cut -f1-6 > ${SAMPLEID}.novoalign.tag.uniq.del.bed; fi
if [ ! -f ${SAMPLEID}.novoalign.tag.uniq.sub.bed ]; then awk '{if($9==">") {print $0}}' ${SAMPLEID}.novoalign.tag.uniq.mutation.txt | cut -f1-6 > ${SAMPLEID}.novoalign.tag.uniq.sub.bed; fi
if [ ! -f ${SAMPLEID}.novoalign.tag.uniq.ins.bed ]; then awk '{if($9=="+") {print $0}}' ${SAMPLEID}.novoalign.tag.uniq.mutation.txt | cut -f1-6 > ${SAMPLEID}.novoalign.tag.uniq.ins.bed; fi

# delte the cache if previously run
rm -rf cache*

# repeat these steps for each of the above..
for i in del sub ins; do
 
  #  Prepare a bedGraph file of mutations in unique tags for visualization in the genome browser (note that the Bedgraph 
  # file of unique CLIP tags is prepared in Step 98). Substitution and insertion files can be converted similarly if necessary.
  if [ ! -f ${SAMPLEID}.novoalign.tag.uniq.${i}.bedgraph ]; then
      perl $CIMSdir/tag2profile.pl -v -ss -exact -of bedgraph -n "Deletion Profile" ${SAMPLEID}.novoalign.tag.uniq.${i}.bed ${SAMPLEID}.novoalign.tag.uniq.${i}.bedgraph
  fi


  # Cluster deletions and evaluate the reproducibility of clustered deletion sites by permutation:
  if [ ! -f ${SAMPLEID}.novoalign.tag.uniq.${i}.CIMS.txt ]; then
    perl $CIMSdir/CIMS.pl -v -n 5 -p -c ./cache_${i} --keep-cache ${SAMPLEID}.novoalign.tag.uniq.bed ${SAMPLEID}.novoalign.tag.uniq.${i}.bed ${SAMPLEID}.novoalign.tag.uniq.${i}.CIMS.txt
    # fi

    # Get the number of ${i}etions in each position relative to the 5′ end of reads (or equivalently the CLIP tag):
    sort -n ./cache_${i}/mutation.pos.txt | uniq -c > num_${i}_mutations.txt
  fi

  # Select robust CIMS with FDR ≤0.001:
  if [ ! -f ${SAMPLEID}.novoalign.tag.uniq.${i}.CIMS.s30.txt ]; then
  echo "Finding robust CIMS"
  awk '{if($9<=0.001) {print $0}}' ${SAMPLEID}.novoalign.tag.uniq.${i}.CIMS.txt | sort -k 9,9n -k 8,8nr -k 7,7n > ${SAMPLEID}.novoalign.tag.uniq.${i}.CIMS.s30.txt
  fi
  # This command will also sort CIMS first by FDR (ascending), then m (descending) and then k (ascending), so that the most reliable sites are ranked at the top, which will facilitate visual inspection of sites in the genome browse


  # Specify −20 to +20 nt around CIMS in a BED file for motif analysis. The exact range of sequences depends on the specific RNABPs, but the ranges from (−10, +10) to (−30, +30) will be a good start.
  if [ ! -f ${SAMPLEID}.novoalign.tag.uniq.${i}.CIMS.s30.41nt.bed ]; then
    echo "Make a bed file for the CIMS +/-20nt for motif analysis"
    awk '{print $1"\t"$2-20"\t"$3+20"\t"$4"\t"$5"\t"$6}' ${SAMPLEID}.novoalign.tag.uniq.${i}.CIMS.s30.txt > ${SAMPLEID}.novoalign.tag.uniq.${i}.CIMS.s30.41nt.bed
  fi
  if [ ! -f ${SAMPLEID}.novoalign.tag.uniq.${i}.CIMS.s30.101nt.bed ]; then
    echo "Make a bed file for the CIMS +/-50nt for motif analysis"
    awk '{print $1"\t"$2-50"\t"$3+50"\t"$4"\t"$5"\t"$6}' ${SAMPLEID}.novoalign.tag.uniq.${i}.CIMS.s30.txt > ${SAMPLEID}.novoalign.tag.uniq.${i}.CIMS.s30.101nt.bed
  fi
  
  # make a fasta of the top CIMs regions for motif analysis
  tmp_41ntbed="${SAMPLEID}.novoalign.tag.uniq.${i}.CIMS.s30.41nt.bed"
  tmp_101ntbed="${SAMPLEID}.novoalign.tag.uniq.${i}.CIMS.s30.101nt.bed"
  tmp_41ntfasta="${tmp_41ntbed%%.bed}.fa"
  tmp_101ntfasta="${tmp_101ntbed%%.bed}.fa"
  echo "Generating fasta files from CIMS bed files."
  module unload bedtools
  module unload gcc
  module load bedtools/2.22.0

  bedtools getfasta -s -fi "$GENOMEFA" -bed "$tmp_41ntbed" -fo "$tmp_41ntfasta"
  bedtools getfasta -s -fi "$GENOMEFA" -bed "$tmp_101ntbed" -fo "$tmp_101ntfasta"

done

# everything after this is optional and should be customized as needed!

# Annotate all peaks files with HOMER; use all BED files, annotate into new dir for each
tmpPWD="$PWD"
tmpGENOME=$(basename $GENOME)
for i in *.bed; do
  # set the outdir
  tmp_anno_dir="${i}_annotation"
  # delete it if it already exists
  rm -rf "$tmp_anno_dir"
  # only annotate certain peak files
  if [[ ! $i == *"101nt"* ]] && [[ ! $i == *"cluster.PH"* ]] && [[ ! $i == *"tag.bed"* ]] && [[ ! $i == *".collapsed.bed"* ]]
  then
    # cd "$tmpPWD"
    echo -e "Annotating file:\n${i}"
    # make the outdir
    mkdir -p "$tmp_anno_dir"
    # annotate
    annotatePeaks.pl "$i" "$tmpGENOME" -annStats ${tmp_anno_dir}/annotation_stats.txt > ${tmp_anno_dir}/annotated_peaks.txt
    # get the number of types of peaks
    Peak_Type_Stats=${tmp_anno_dir}/peak_type_stats.txt
    tail -n +2 ${tmp_anno_dir}/annotated_peaks.txt | cut -f8 | cut -d '(' -f1 | sort | uniq -c | sed -e 's/ *//' -e 's/ /\t/' -e "s/'//" -e 's/NA/Other/' -e 's/ //' > $Peak_Type_Stats
  fi
done





# exit

# MEME-ChIP motif analysis
module unload gcc
module unload perl
module unload python
module load perl/5.22.1

module unload meme/4.9.1
module load meme/4.11.1
tmp_JasparDB="/ifs/home/kellys04/data/motifs/motif_databases/JASPAR/JASPAR_CORE_2016_vertebrates.meme"
tmp_JOLMAdb="/ifs/home/kellys04/data/motifs/motif_databases/EUKARYOTE/jolma2013.meme"
# tmpMEME_DB2, tmpMEME_DB
# remember the present working directory..
tmpPWD="$PWD"
for i in del sub ins; do
  for q in 41 101; do
    # for k in 1 2 3; do # dont need to repeat this several times, MEME gives consistent output!
      # repeat 3 times because with low # input seq's you get some random variance
      # reset the PWD
      cd "$tmpPWD"
      
      # tmp_motifdir="motifsMEME-CHIP-2/${i}_${q}nt_rep${k}"
      tmp_motifdir="motifsMEME-CHIP-2/${i}_${q}nt"
      if [ ! -d $tmp_motifdir ]; then 
        echo -e "Motif outdir is\n ${tmp_motifdir}"
        mkdir -p "$tmp_motifdir"
        tmpGENOME=$(basename $GENOME)
        input_bed="${SAMPLEID}.novoalign.tag.uniq.${i}.CIMS.s30.${q}nt.bed"
        tmp_fasta="${input_bed%%.bed}.fa"
        # "$tmpMEME_DB"
        # -db "$$tmpMEME_DB"
        meme-chip -oc "$tmp_motifdir" -index-name meme-chip.html -time 300 -order 1 -norand -db "$tmp_JOLMAdb" -db "$tmp_JasparDB" -db "$tmpMEME_DB2" -db "$tmpMEME_DB" -meme-mod zoops -meme-minw 6 -meme-maxw 30 -meme-nmotifs 5 -dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 "$tmp_fasta"
        # meme-chip -oc "$tmp_motifdir" -index-name meme-chip.html -time 300 -meme-p 4 -meme-nmotifs 5  -db "$tmp_JasparDB" "$tmp_fasta"
      fi
    # done
  done
done


# HOMER motif analysis
# remember the present working directory..
tmpPWD="$PWD"

for i in del sub ins; do
  for q in 41 101; do
    for k in 1 2 3; do
      # repeat 3 times because with low # input seq's you get some random variance
      # reset the PWD
      cd "$tmpPWD"

      # try motif analysis with HOMER
      # make an outdir, delete if it already exists & switch to it
      tmp_motifdir="motifsHOMER_collapse/${i}_${q}nt_rep${k}"
      if [ ! -d $tmp_motifdir ]; then 
        echo -e "Motif outdir is\n ${tmp_motifdir}"
        # rm -rf "$tmp_motifdir"
        mkdir -p "$tmp_motifdir"
        # cd "$tmp_motifdir"

        tmpGENOME=$(basename $GENOME)
        # PREPARSED_DIR="$HOME/software/homer/preparsed/${tmpGENOME}"
        PREPARSED_DIR="$tmp_motifdir/preparsed"
        # -preparse
        # 
        
        input_bed="${SAMPLEID}.novoalign.tag.uniq.${i}.CIMS.s30.${q}nt.bed"
        collapse_bed="${SAMPLEID}.novoalign.tag.uniq.${i}.CIMS.s30.${q}nt.collapsed.bed"
        # need to collapse the CIMS for HOMER motif analysis
        # mergePeaks peakFile.txt -strand > output.peaks.txt
        if [ ! -f $collapse_bed ]; then mergePeaks "$input_bed" -strand > "$collapse_bed"; fi
        
        # run HOMER motif analysis
        # echo "Running HOMER for motif analysis"
        start_time1=`date +%s`
        # findMotifsGenome.pl ${SAMPLEID}.novoalign.tag.uniq.${i}.CIMS.s30.${q}nt.bed $tmpGENOME "$tmp_motifdir" -preparsedDir $PREPARSED_DIR -preparse -size given -p $THREADS -dumpFasta -rna
        findMotifsGenome.pl "$collapse_bed" $tmpGENOME "$tmp_motifdir" -preparsedDir $PREPARSED_DIR -preparse -size given -p $THREADS -dumpFasta -rna
        echo -e "Execution time:\t$(expr $(date +%s) - $start_time1)s"
        # start_time1=`date +%s`


        # start_time1=`date +%s`
        # annotatePeaks.pl "${SAMPLEID}.novoalign.tag.uniq.${i}.CIMS.s30.${q}nt.bed" "$tmpGENOME" -annStats $tmp_motifdir/annotation_stats.txt > $tmp_motifdir/annotated_peaks.txt
        annotatePeaks.pl "$collapse_bed" "$tmpGENOME" -annStats $tmp_motifdir/annotation_stats.txt > $tmp_motifdir/annotated_peaks.txt
        # echo -e "Execution time:\t$(expr $(date +%s) - $start_time1)s"
        #-genomeOntology $ANNOT_OUT_DIR/genomeOntology
        # -d $ANNOT_OUT_DIR/tag_directory

        # get the number of types of peaks
        Peak_Type_Stats=$tmp_motifdir/peak_type_stats.txt

        # get the unique entries from the types of peak annotations, output in tab delim format
        tail -n +2 $tmp_motifdir/annotated_peaks.txt | cut -f8 | cut -d '(' -f1 | sort | uniq -c | sed -e 's/ *//' -e 's/ /\t/' -e "s/'//" -e 's/NA/Other/' -e 's/ //' > $Peak_Type_Stats
      fi
  done
  done
done
