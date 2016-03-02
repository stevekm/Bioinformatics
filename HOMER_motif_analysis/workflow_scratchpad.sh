#!/bin/bash

# this is where I set, edit, and tweak parameters in order to fit the exact project I am working on (e.g. file structure, filenames, etc.)
# use this as a scratch-pad to set up work, then copy/paste into the bash shell to run the commands

# change to proj dir
ProjDir="/ifs/home/$(whoami)/projects/project_directory"
cd $ProjDir
Results_Share_Dir="/ifs/home/$(whoami)/projects/project_directory/results_share_dir"
# dir that contains nothing but ".bed" and "_scores.bed" files for each sample
Results_Peaks_Dir="$Results_Share_Dir/peaks"

# where to put the output of the analysis
Motif_Analysis_Outdir="$ProjDir/motif_analysis"

# primary scripts for qsub submission with commands to 
Motif_Analysis_Script="/ifs/home/$(whoami)/projects/project_directory/code/motif_analysis_HOMER3.sh"
Peak_Annotation_Script="/ifs/home/$(whoami)/projects/project_directory/code/peak_annotation_HOMER.sh"
Tag_Dir_Script="/ifs/home/$(whoami)/projects/project_directory/code/make_tagDirectory_HOMER-bed.sh"
# make sure the scripts are executable 
chmod +x "$Motif_Analysis_Script"
chmod +x "$Peak_Annotation_Script"
chmod +x "$Tag_Dir_Script"

# commonly used HOMER parameters go here
Motif_Analysis_Script_params="/ifs/home/$(whoami)/projects/project_directory/code/motif_analysis_HOMER2_params.sh"


#~~~ Shell Options ~~~~#
# Steve's favorite glob settings
# remember whether extglob was originally set, so we know whether to unset it
shopt -q extglob; extglob_set=$?
# set extglob if it wasn't originally set.
((extglob_set)) && shopt -s extglob
# Note, 0 (true) from shopt -q is "false" in a math context.
shopt -q nullglob; nullglob_set=$?
((nullglob_set)) && shopt -s nullglob
shopt -q globstar; globstar_set=$?
((globstar_set)) && shopt -s globstar
#~~~~~~~~~~~~~~~~~~#


#~~~ Sample Workflow Code ~~~~#
# this code is changed to match the project being worked on; your basic code goes here
# figure out what you need to do to get the shell to grab each .bed file and pass it to the analysis script

# glob up the bed files for the peaks
# problem: We want to get all the files that end with '.bed' but not the files that end with '_scores.bed'
# Results_Share_Dir_Peaks=($Results_Share_Dir/peaks/**/!(*_scores.bed)) # this almost works but includes parent dirs..
Results_Share_Dir_BEDs=($Results_Peaks_Dir/**/@(*[^_scores].bed))

for i in "${Results_Share_Dir_BEDs[@]}"; do
  echo "$i"
  
  # return the base file and its parent dir only, by stripping the dir path substring pattern from it
  Dir_File_base="${i#$Results_Peaks_Dir/}"
  
  # type = name of .bed file parent dir
  tmp_Type=$(echo $Dir_File_base | cut -d '/' -f1)
  
  # file = .bed file
  tmp_File=$(basename "$i")
  
  # sample = filename minus file extension
  tmp_SampleID=${tmp_File%.*}
  
  # output dir for the sample's motif analysis
  tmp_motif_outdir="${Motif_Analysis_Outdir}"/"${tmp_Type}"/"${tmp_SampleID}"
  
  # verbose output for testing
  echo "$tmp_Type is tmp_Type"
  echo "$tmp_File is tmp_File"
  echo "$tmp_SampleID is tmp_SampleID"
  echo "$tmp_motif_outdir is tmp_motif_outdir"
  
  # make the outdir
  mkdir -p "$tmp_motif_outdir"
  
  # start the completion timer # I didn't actually test this timer it might not work yet
  Job_Timer_File="${tmp_motif_outdir}"/job_timer.${RANDOM}.txt
  echo -e "job start: $(date +[%F]%T)\n" > $Job_Timer_File
  echo -e "command is\n" >> $Job_Timer_File
  echo -e "qsub -q all.q -wd $tmp_motif_outdir -o :${tmp_motif_outdir}/ -e :${tmp_motif_outdir}/ -pe threaded 8-32 $Motif_Analysis_Script ${tmp_motif_outdir} ${i} ${Motif_Analysis_Script_params}" >> $Job_Timer_File
  
  
  
  # submit the analysis job to qsub
  # make Tag Directory
  qsub -q all.q -wd $tmp_TagDir_outdir -o :${tmp_TagDir_outdir}/ -e :${tmp_TagDir_outdir}/ -pe threaded 4 "$Tag_Dir_Script" "${tmp_TagDir_outdir}" "${i}" "${Motif_Analysis_Script_params}"
  
  # Annotate Peaks
  qsub -q all.q -wd $tmp_PeakAnnot_outdir -o :${tmp_PeakAnnot_outdir}/ -e :${tmp_PeakAnnot_outdir}/ -pe threaded 4 "$Peak_Annotation_Script" "${tmp_PeakAnnot_outdir}" "${i}" "${Motif_Analysis_Script_params}"
  
  # Motif Analysis
 qsub -q all.q -wd $tmp_motif_outdir -o :${tmp_motif_outdir}/ -e :${tmp_motif_outdir}/ -pe threaded 8-24 "$Motif_Analysis_Script" "${tmp_motif_outdir}" "${i}" "${Motif_Analysis_Script_params}"
  # # for qsub; 
  # -q all.q : selects only nodes with name "all.q*"; excludes GPU, HighMem nodes
  # -wd $tmp_motif_outdir : sets the working directory for the script to execute from
  # -o :${tmp_motif_outdir}/ : sets the outdir for std out logs
  # -e :${tmp_motif_outdir}/ : sets the outdir for the std err logs
  # -pe threaded 8-32 : sets the range of threads allotted to the job
  
  
  # complete the timer
  echo -e "job finish: $(date +[%F]%T)\n" >> $Job_Timer_File
done
# change back to proj dir
# cd $ProjDir

#~~~ Shell Options ~~~~#
# unset globs if it wasn't originally set
((extglob_set)) && shopt -u extglob
((nullglob_set)) && shopt -u nullglob
((globstar_set)) && shopt -u globstar
#~~~~~~~~~~~~~~~~~~#

