#!/bin/bash

# this workflow converts Diffbind output csv files into GREAT format BED files

GREAT_dir="$HOME/projects/SmithLab_2016-03-10/project_notes/GREAT"
tmp_diffbind_dir="/ifs/data/sequence/results/smithlab/2016-01-04/ChIP-Seq/diffbind"
cd $GREAT_dir

# # Diffbind subdirs per histone mark
# ll $tmp_diffbind_dir
# total 206K
# drwxrws--- 7 kellys04 bin 125 Mar 18 15:16 .
# drwxrws--- 8 kellys04 bin 219 Mar 18 19:49 ..
# drwxrws--- 2 id460    bin 359 Mar 18 15:12 H3K27AC
# drwxrws--- 2 id460    bin 307 Mar 18 14:46 H3K27ME3
# drwxrws--- 2 id460    bin 402 Mar 18 14:55 H3K4ME3
# drwxrws--- 2 id460    bin 316 Mar 18 14:32 H3K9AC
# drwxrws--- 2 id460    bin 212 Mar 18 14:41 H3K9ME3

# # sample Diffbind files
# ll ${tmp_diffbind_dir}/H3K27AC
# drwxrws--- 2 id460    bin  359 Mar 18 15:12 .
# drwxrws--- 7 kellys04 bin  125 Mar 18 15:16 ..
# -rw-rw---- 1 id460    bin  92K Mar 18 15:12 diff_bind.D-vs-R.blocking.p005.csv
# -rw-rw---- 1 id460    bin 944K Mar 18 15:12 diff_bind.D-vs-R.blocking.p020.csv
# -rw-rw---- 1 id460    bin  77M Mar 18 15:12 diff_bind.D-vs-R.blocking.p100.csv
# -rw-rw---- 1 id460    bin  23K Mar 18 14:44 diff_bind.D-vs-R.p020.csv
# -rw-rw---- 1 id460    bin  77M Mar 18 14:44 diff_bind.D-vs-R.p100.csv
# -rw-rw---- 1 id460    bin 6.7K Mar 18 13:57 diffbind-sample-sheet.csv
# -rw-rw---- 1 id460    bin  17K Mar 18 14:33 plot.heatmap.fpkm.pdf
# -rw-rw---- 1 id460    bin 7.6K Mar 18 14:33 plot.pca.fpkm.pdf

# # contents of a Diffbind file
# head -n 2 ${tmp_diffbind_dir}/H3K27AC/diff_bind.D-vs-R.blocking.p005.csv
# "seqnames","start","end","width","strand","Conc","Conc_D","Conc_R","Fold","p.value","FDR","AGK.D.H3K27AC","BVI.D.H3K27AC","CBK.D.H3K27AC","DKJ.D.H3K27AC","FLV.D.H3K27AC","GHW.D.H3K27AC","IDY.D.H3K27AC","PRW.D.H3K27AC","SPN.D.H3K27AC","SRR.D.H3K27AC","ZGR.D.H3K27AC","ZNK.D.H3K27AC","AGK.R.H3K27AC","BVI.R.H3K27AC","CBK.R.H3K27AC","DKJ.R.H3K27AC","FLV.R.H3K27AC","GHW.R.H3K27AC","IDY.R.H3K27AC","PRW.R.H3K27AC","SPN.R.H3K27AC","SRR.R.H3K27AC","ZGR.R.H3K27AC","ZNK.R.H3K27AC","feature","external_gene_name","gene_biotype","start_position","end_position","insideFeature","distancetoFeature","shortestDistance","fromOverlappingOrNearest"
# "chr1",3229502,3234925,5424,"*",5.32490538191405,4.51171079235395,5.84180781927711,-1.33009702692316,2.99943604239429e-05,0.0367619699255697,89.3284837621687,1.08112603875298,0.782677655086077,123.602132308718,3.26195148006001,0.756951969765614,6.86268034787737,22.6649228252718,16.8275215043467,0.624461805628116,3.52220337394804,4.42694388057497,184,18.1914591885298,3.00040740592174,101.040978906413,31.2572176942067,49.105436537724,32.8178520162554,66.3389402951096,22.5330310513055,30.2513039459725,5.39137251101084,144.313139334195,"ENSG00000142611","PRDM16","protein_coding",2985732,3355185,"inside",246482,120260,"shortestDistance"

# tail -n +2 diff_bind.D-vs-R.blocking.p100.csv | head | cut -d ',' -f1-4

# just get the coordinates of the file, save BED format
mkdir -p ${GREAT_dir}/${tmp_mark}
tail -n +2 ${tmp_diffbind_dir}/${tmp_mark}/${tmp_file} | cut -d ',' -f1-3 | sed -e 's/,/\t/g' -e 's/"//g' > ${GREAT_dir}/${tmp_mark}/${tmp_file}


for i in H3K27AC H3K27ME3 H3K4ME3 H3K9AC H3K9ME3; do
  # get the CSV files to be converted
  FILES=${tmp_diffbind_dir}/${i}/*.csv
  echo "$i"
  
  # make the outdir
  mkdir -p  ${GREAT_dir}/${i}
  
  # convert each file
  for q in $FILES; do
    # get the filename, convert to bed file
    tmp_filename=$(basename $q)
    tmp_filename=${tmp_filename/.csv/.bed}
    echo "$tmp_filename"
    
    # strip headers, convert , to tab delimiter, remove the quotes, only first 3 columns
    tail -n +2 "$q" | cut -d ',' -f1-3 | sed -e 's/,/\t/g' -e 's/"//g' > ${GREAT_dir}/${i}/${tmp_filename}
    done
done

# # contents of files after processing
# head -n 5 ${GREAT_dir}/H3K27AC/diff_bind.D-vs-R.blocking.p005.bed
# chr1	3229502	3234925
# chr1	6974643	6975743
# chr1	6976806	6979695
# chr1	9335331	9336936
# chr1	17685135	17687663
