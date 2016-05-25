#!/bin/bash

# a workflow script that performs HOMER mergePeaks on 4 sets of peaks (2 TSS regions, 2 histone marks) to find the peaks in common

# Ensembl TSS region bed file
ens_bed="$HOME/projects/SmithLab_ChIPSeq_2016-12-31/project_data/genes.ensembl.GRCh37.82_noGLMT_chr.tss_collapse_10kb.bed"
# Gencode TSS region bed file
gen_bed="$HOME/projects/SmithLab_ChIPSeq_2016-12-31/project_data/gencode.v19.annotation.tss_collapse_10kb.bed"

# iterate over list of types of peaks to use; these are components of the filepath
for q in sample group; do
# iterator; type of peak
tmp_type="$q"

# load the correct version of HOMER
module load homer/v4.6

# project directory
ProjDir="$HOME/projects/SmithLab_ChIPSeq_2016-12-31"
# path to parent dir containing peak calling results
peaks_resultsDir="${ProjDir}/pipeline/peaks/results"
# outdir for peak overlaps
testDir="${ProjDir}/project_notes/peak_overlap/macs_broad_by_${tmp_type}_10kbp_direct2"
# project sample sheet
samplesheet="${ProjDir}/inputs/sample-sheet.tsv"

# make sure the outdir exists, cd to it
mkdir -p "$testDir"
cd "$testDir"

# iterate over a list of the SampleID; get this from the 9th field in the Samplesheet!
# # example SampleID: ABC, EGF, FAS
for i in $(tail -n +2 $samplesheet | cut -f9 | sort -u); do
# ( # use this for running them all at once, be careful with this! See matching ) & done at the end of the loop

  # get the sample ID, setup the dir
  echo  "$i"
  tmp_sampleID="$i"
  
  # make a subdir for each sample output, cd to it
  tmp_outDir="${testDir}/${tmp_sampleID}"
  mkdir -p "$tmp_outDir"
  cd "$tmp_outDir"
  
  
  
  
  # find the peak files for each histone mark, copy to the new outdir
  # # set output files to copy to
  tmp_bedfile_H3K27AC="${tmp_outDir}/${tmp_sampleID}_R_H3K27AC.bed"
  echo "tmp_bedfile_H3K27AC is $tmp_bedfile_H3K27AC"
  
  tmp_bedfile_H3K4ME3="${tmp_outDir}/${tmp_sampleID}_R_H3K4ME3.bed"
  echo "tmp_bedfile_H3K4ME3 is $tmp_bedfile_H3K4ME3"
  
  find "$peaks_resultsDir" -type f -name "peaks.bed" -path "*/peaks.by_${tmp_type}.macs_broad/*" -path "*/${tmp_sampleID}-R-H3K27AC*" -exec cp {} "${tmp_bedfile_H3K27AC}" \;
  
  find "$peaks_resultsDir" -type f -name "peaks.bed" -path "*/peaks.by_${tmp_type}.macs_broad/*" -path "*/${tmp_sampleID}-R-H3K4ME3*" -exec cp {} "${tmp_bedfile_H3K4ME3}" \;
  
  
  # make sure both histone marks had peaks before trying to continue
  if [[ -f ${tmp_bedfile_H3K4ME3} && -f ${tmp_bedfile_H3K27AC} ]]; then
    echo "Both files found"
    
    # copy the reference TSS regions; makes overlapping easier
    cp "$ens_bed" "ensemble.bed"
    cp "$gen_bed" "gencode.bed"
    
    tmp_ens_bed="ensemble.bed"
    tmp_gen_bed="gencode.bed"
    tmp_outH3K4ME3="$(basename ${tmp_bedfile_H3K4ME3})"
    tmp_outH3K27AC="$(basename ${tmp_bedfile_H3K27AC})"
    
    # use HOMER mergePeaks to overlap the peaks in the bed files
    # # prefix and file order need to stay consistent throughout, for the rest of the script to work
    mergePeaks "$tmp_outH3K4ME3" "$tmp_outH3K27AC" "$tmp_ens_bed" "$tmp_gen_bed" -prefix mergepeaks -venn venn.txt -matrix matrix.txt
    
    # print console output
    echo -e "\n"
    cat matrix.txt.count.matrix.txt
    echo -e "\n"
    cat venn.txt
    echo -e "\n"
    
    # set the mergePeaks outputs names
    tmp_mergeH3K4ME3="mergepeaks_${tmp_outH3K4ME3}"
    tmp_mergeH3K27AC="mergepeaks_${tmp_outH3K27AC}"
      
    # count the unique peaks
    num_H3K4ME3=$(tail -n +2 $tmp_mergeH3K4ME3 | wc -l)
    echo "num_H3K4ME3 is $num_H3K4ME3"
    num_H3K27AC=$(tail -n +2 $tmp_mergeH3K27AC | wc -l)
    echo "num_H3K27AC is $num_H3K27AC"
    num_ensemble=$(tail -n +2 mergepeaks_${tmp_ens_bed} | wc -l )
    echo "num_ensemble is $num_ensemble"
    num_gencode=$(tail -n +2 mergepeaks_${tmp_gen_bed} | wc -l )
    echo "num_gencode is $num_gencode"
    
    # count the peaks in common
    num_H3K27AC_ens=$(tail -n +2 "mergepeaks_${tmp_outH3K27AC}_${tmp_ens_bed}" | wc -l)
    echo "num_H3K27AC_ens is $num_H3K27AC_ens"
    num_H3K27AC_gen=$(tail -n +2 "mergepeaks_${tmp_outH3K27AC}_${tmp_gen_bed}" | wc -l)
    echo "num_H3K27AC_gen is $num_H3K27AC_gen"
    num_H3K4ME3_ens=$(tail -n +2 "mergepeaks_${tmp_outH3K4ME3}_${tmp_ens_bed}" | wc -l)
    echo "num_H3K4ME3_ens is $num_H3K4ME3_ens"
    num_H3K4ME3_gen=$(tail -n +2 "mergepeaks_${tmp_outH3K4ME3}_${tmp_gen_bed}" | wc -l)
    echo "num_H3K4ME3_gen is $num_H3K4ME3_gen"
    
    num_H3K4ME3_H3K27AC=$(tail -n +2 "mergepeaks_${tmp_outH3K4ME3}_${tmp_outH3K27AC}" | wc -l)
    echo "num_H3K4ME3_H3K27AC is $num_H3K4ME3_H3K27AC"
    num_ens_gen=$(tail -n +2 "mergepeaks_${tmp_ens_bed}_${tmp_gen_bed}" | wc -l)
    echo "num_ens_gen is $num_ens_gen"
    
    num_H3K27AC_ens_gen=$(tail -n +2 "mergepeaks_${tmp_outH3K27AC}_${tmp_ens_bed}_${tmp_gen_bed}" | wc -l)
    echo "num_H3K27AC_ens_gen is $num_H3K27AC_ens_gen"
    num_H3K4ME3_ens_gen=$(tail -n +2 "mergepeaks_${tmp_outH3K4ME3}_${tmp_ens_bed}_${tmp_gen_bed}" | wc -l)
    echo "num_H3K4ME3_ens_gen is $num_H3K4ME3_ens_gen"
    
    
    
    num_H3K4ME3_H3K27AC_ens=$(tail -n +2 "mergepeaks_${tmp_outH3K4ME3}_${tmp_outH3K27AC}_${tmp_ens_bed}" | wc -l)
    echo "num_H3K4ME3_H3K27AC_ens is $num_H3K4ME3_H3K27AC_ens"
    
    num_H3K4ME3_H3K27AC_gen=$(tail -n +2 "mergepeaks_${tmp_outH3K4ME3}_${tmp_outH3K27AC}_${tmp_gen_bed}" | wc -l)
    echo "num_H3K4ME3_H3K27AC_ens is $num_H3K4ME3_H3K27AC_ens"

    
    num_H3K4ME3_H3K27AC_ens_gen=$(tail -n +2 "mergepeaks_${tmp_outH3K4ME3}_${tmp_outH3K27AC}_${tmp_ens_bed}_${tmp_gen_bed}" | wc -l)
    echo "num_H3K4ME3_H3K27AC_ens_gen is $num_H3K4ME3_H3K27AC_ens_gen"
    
    # NOTE: If a comparison yielded no peaks, no file will be output; set this as 0 peaks for Venn when passing values to R for plotting    
    
    # load correct version of R; need 3.2.0 for compatibility with VennDiagram package as of April 26, 2016  
    module unload r
    module load r/3.2.0
    # pass values to R for plotting; if doesn't exist, pass 0 instead
    # # use Rscript heredoc so we don't need a separate R script file!
    Rscript --slave --no-save --no-restore - "$tmp_sampleID" "${num_H3K4ME3:=0}" "${num_H3K27AC:=0}" "${num_ensemble:=0}" "${num_gencode:=0}" "${num_H3K27AC_ens:=0}" "${num_H3K27AC_gen:=0}" "${num_H3K4ME3_ens:=0}" "${num_H3K4ME3_gen:=0}" "${num_H3K4ME3_H3K27AC:=0}" "${num_ens_gen:=0}" "${num_H3K27AC_ens_gen:=0}" "${num_H3K4ME3_ens_gen:=0}" "${num_H3K4ME3_H3K27AC_ens:=0}" "${num_H3K4ME3_H3K27AC_gen:=0}" "${num_H3K4ME3_H3K27AC_ens_gen:=0}" <<EOF
      ## R code
      # load packages 
      library('VennDiagram')
      library('gridExtra')
      # print to console
      cat("\nR loaded\n")
      # get the script arguments
      args <- commandArgs(TRUE); cat("Script args are:\n"); args
      SampleID<-args[1]
      peaks_H3K4ME3<-as.numeric(args[2]) # num_H3K4ME3
      peaks_H3K27AC<-as.numeric(args[3]) # num_H3K27AC
      peaks_ensembl<-as.numeric(args[4]) # num_ensemble
      peaks_gencode<-as.numeric(args[5]) # num_gencode
      peaks_H3K27AC_ens<-as.numeric(args[6]) # num_H3K27AC_ens
      peaks_H3K27AC_gen<-as.numeric(args[7]) # num_H3K27AC_gen
      peaks_H3K4ME3_ens<-as.numeric(args[8]) # num_H3K4ME3_ens
      peaks_H3K4ME3_gen<-as.numeric(args[9]) # num_H3K4ME3_gen
      peaks_H3K4ME3_H3K27AC<-as.numeric(args[10]) # num_H3K4ME3_H3K27AC
      peaks_ens_gen<-as.numeric(args[11]) # num_ens_gen
      peaks_H3K27AC_ens_gen<-as.numeric(args[12]) # num_H3K27AC_ens_gen
      peaks_H3K4ME3_ens_gen<-as.numeric(args[13]) # num_H3K4ME3_ens_gen
      peaks_H3K4ME3_H3K27AC_ens<-as.numeric(args[14]) # num_H3K4ME3_H3K27AC_ens
      peaks_H3K4ME3_H3K27AC_gen<-as.numeric(args[15]) # num_H3K4ME3_H3K27AC_gen
      peaks_H3K4ME3_H3K27AC_ens_gen<-as.numeric(args[16]) # num_H3K4ME3_H3K27AC_ens_gen
      
      # convert to venn ranges
      # this looks redundant but its easier to understand this way!
      n1<-peaks_H3K4ME3
      n2<-peaks_H3K27AC
      n3<-peaks_ensembl
      n4<-peaks_gencode
      
      n12<-peaks_H3K4ME3_H3K27AC
      n13<-peaks_H3K4ME3_ens
      n14<-peaks_H3K4ME3_gen
      
      n23<-peaks_H3K27AC_ens
      n24<-peaks_H3K27AC_gen
      
      n34<-peaks_ens_gen
      
      n134<-peaks_H3K4ME3_ens_gen
      n234<-peaks_H3K27AC_ens_gen
      
      n123<-peaks_H3K4ME3_H3K27AC_ens
      n124<-peaks_H3K4ME3_H3K27AC_gen
      
      n1234<-peaks_H3K4ME3_H3K27AC_ens_gen
      
      # set filename for plot, based on SampleID
      plot_filename<-paste0(SampleID,"_10kb_direct_TSS_peaks.pdf") 


      # set up the 4-way Venn object, don't print it yet
      venn<-draw.quad.venn(
        area1=n1+n1234+n123+n124+n12+n134+n13+n14,
        area2=n2+n1234+n123+n124+n12+n234+n23+n24,
        area3=n3+n1234+n123+n134+n13+n234+n23+n34,
        area4=n4+n1234+n124+n134+n14+n234+n24+n34,
        n12=n12+n1234+n123+n124,
        n13=n13+n1234+n123+n134,
        n14=n14+n1234+n124+n134,
        n23=n23+n1234+n123+n234,
        n24=n24+n1234+n124+n234,
        n34=n34+n1234+n134+n234,
        n123=n123+n1234,
        n124=n124+n1234,
        n134=n134+n1234,
        n234=n234+n1234,
        n1234=n1234,
        category=c('H3K4ME3','H3K27AC',"Ensembl","GenCode"),
        fill=c("red","blue","yellow","green"),
        alpha=c(0.3,0.3,0.3,0.3),
        cex=2,
        cat.cex=2,
        ind=FALSE)
      
      # print the plot to a PDF, include a title at the top
      pdf(plot_filename,width = 8,height = 8)
      grid.arrange(gTree(children=venn), top=paste0(SampleID," TSS 10kbp Region Direct Overlap")) #, bottom="subtitle")
      dev.off()
EOF
    else
    echo "A file is missing!"
  fi
  echo -e "\n\n"
  cd "$testDir"
#) & done # # use this for running them all at once, be careful with this!
 done

# for simplicity, copy all the PDFs to a common directory as well; this won't work the first time if you ran all of them at once
mkdir -p "${testDir}/all_pdf_10kbp_direct"
find . -type f -name "*.pdf" -exec cp {} "${testDir}/all_pdf_10kbp_direct/" \;

done
