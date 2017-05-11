#!/bin/bash

# Workflow for running ROSE Super Enhancer calling on HiC-Bench pipeline output

# source:
#
# git clone https://stevekm@bitbucket.org/young_computation/rose.git
# https://bitbucket.org/young_computation/rose/overview
# http://younglab.wi.mit.edu/super_enhancer_code.html

# - only use H3K27AC peaks, MACS broad, per sample

# ~~~~~~~~~~ SETTINGS & LOCATIONS ~~~~~~~~~~ #
se_dir="/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2016-12-31/project_notes/SuperEnhancer_calling"
se_outdir="${se_dir}/output"; mkdir -p "$se_outdir"
samplesheet="/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2016-12-31/inputs/sample-sheet.tsv"
pipeline_dir="/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2016-12-31/pipeline"
rose_dir="/ifs/home/kellys04/software/rose" # place where I cloned the ROSE repo
div="--------------------"


# ~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~ #
find_bam () {
    local sampleID="$1"
    local align_dir="${pipeline_dir}/align/results"
    # pipeline/align/results/align.by_sample.bowtie2/control-3pos/alignments.bam
    find "$align_dir" -path "*/${sampleID}/*" -name "alignments.bam"
}

find_bai () {
    local sampleID="$1"
    local align_dir="${pipeline_dir}/align/results"
    # pipeline/align/results/align.by_sample.bowtie2/K27Ac-2/alignments.bam.bai
    find "$align_dir" -path "*/${sampleID}/*" -name "alignments.bam.bai"
}

find_bed () {
    local sampleID="$1"
    local peaks_dir="${pipeline_dir}/peaks/results/peaks.by_sample.macs_broad"
    find "$peaks_dir" -path "*/${sampleID}/*" -name "peaks.bed"
}

clean_gff () {
    # remove header and 'bed2gff' label
    # be careful with this step because one time I messed it up and ROSE had a massive memory leak and knocked out several nodes on the cluster
    local input_gff="$1"
    local tmp_file="$(dirname "$input_gff")/tmp"
    cat "$input_gff" | sed 1,3d | sed -e 's/bed2gff[[:space:]]/\t/g' > "$tmp_file" && /bin/mv "$tmp_file" "$input_gff"
}

convert_bed_to_gff () {
    local input_bed_file="$1"
    local output_gff_file="$2"
    local convert_script="/ifs/home/kellys04/software/bed_to_gff_converter.py"
    python ${convert_script} "$input_bed_file" "$output_gff_file"
    clean_gff "$output_gff_file"
}


run_rose () {
    # Need to set up ROSE dir per sample:
    # all sample data files + all ROSE scripts symlinked together in one dir
    local sampleID="$1"
    local control_sample="$2"
    local sample_outdir="${se_outdir}/${sampleID}"
    local rose_outdir="${sample_outdir}/SuperEnhancers"; mkdir -p "$rose_outdir"
    local log_dir="${sample_outdir}/logs"; mkdir -p "$log_dir"
    
    # files needed for ROSE
    local bam_file="$(find_bam "$sampleID")"
    local bai_file="$(find_bai "$sampleID")"
    local bed_file="$(find_bed "$sampleID")"
    
    # symlink filenames to use
    local bam_link_name="${sampleID}_$(basename "$bam_file")"
    local bai_link_name="${sampleID}_$(basename "$bai_file")"
    local bed_link_name="${sampleID}_$(basename "$bed_file")"
    
    # name for the GFF file
    local gff_file="${bed_link_name%%.bed}.gff"
    
    # default argument for ROSE; assume no control supplied
    local ROSE_control_arg=""
    
    # setup sample files
    (
    cd "$sample_outdir"
    ln -fs ${rose_dir}/* ./
    ln -fs "$bam_file" "$bam_link_name"
    ln -fs "$bai_file" "$bai_link_name"
    ln -fs "$bed_file" "$bed_link_name"
    
    # need to convert the BED to GFF format
    convert_bed_to_gff "$bed_link_name" "$gff_file"
    )
    
    # check for a control sample; default value is NA for no control
    if [ "$control_sample" != "NA" ]; then
        echo "Finding control sample files..."
        local control_bam_file="$(find_bam "$control_sample")"
        local control_bai_file="$(find_bai "$control_sample")"
        local control_bam_linkname="${control_sample}_$(basename "$control_bam_file")"
        local control_bai_linkname="${control_sample}_$(basename "$control_bai_file")"
        local ROSE_control_arg="-c $control_bam_linkname"
        # setup control files
        (
        cd "$sample_outdir"
        ln -fs "$control_bam_file" "$control_bam_linkname"
        ln -fs "$control_bai_file" "$control_bai_linkname"
        )
    fi
    
    printf "%s\n%s\n%s\n%s\n%s\n%s\n\n" "$bam_file" "$bai_file" "$bed_file" "$sample_outdir" "$control_bam_file" "$control_bai_file"
    printf "%s\n%s\n%s\n%s\n\n\n" "$bam_link_name" "$bai_link_name" "$bed_link_name" "$gff_file"
    (
    cd "$sample_outdir"
    ls -l
    pwd
    printf "\nLog dir will be:\n%s\n\n" "$log_dir"
    # submit qsub cluster job to run ROSE
     qsub -wd $PWD -o :${log_dir}/ -e :${log_dir}/ -j y -N "$sampleID" <<E0F
set -x
python ./ROSE_main.py -g HG19 -i "$gff_file" -r "$bam_link_name" $ROSE_control_arg -o "$rose_outdir" -t 2500
E0F
    )
}

# ~~~~~~~~~~ RUN ~~~~~~~~~~ #
tail -n +2 "$samplesheet" | while read line; do
    if [ ! -z "$line" ]; then # no empty lines
        sampleID="$(echo "$line" | cut -f1)" 
        case "$sampleID" in 
            *K27Ac*) # sample ID includes K27Ac
            printf "\n%s\n\n" "$div"
            control_sample="$(echo "$line" | cut -f2)"
            echo "$sampleID $control_sample"
            run_rose "$sampleID" "$control_sample"
            ;;
        esac
    fi
done
