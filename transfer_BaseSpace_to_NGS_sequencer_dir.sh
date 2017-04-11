#!/bin/bash

## USAGE: /ifs/data/scripts/transfer_BaseSpace_to_NGS_sequencer_dir.sh "ProjectID" "/path/to/BaseSpace_dir"
## DESCRIPTION: This script will copy .fastq.gz files from a BaseSpace project
## directory to the 'NGS_sequencer_dir' directory. 

# ~~~~~ CHECK SCRIPT ARGS ~~~~~ #
if (( "$#" != "2" )); then
    echo "ERROR: Wrong number of arguments supplied"
    grep '^##' $0
    exit
fi



# ~~~~~ CUSTOM FUNCTIONS ~~~~~ #
copy_to_common_dir () {
    # copy all files to a common output dir
    local source_dir="$1"
    local destination_dir="$2"
    local file_type="$3"
    find "${source_dir}" -name "${file_type}" -exec rsync -vctrPh {} "$destination_dir" \; # --dry-run
}

copy_with_parents () {
    # copy files from BaseSpace BaseMount, keep parent dirs
    # e.g.: /BaseSpace/Projects/NextSeq-123/Samples/ABC-0091/Files/ABC-0091_S6_L002_R2_001.fastq.gz
    # output_path=ABC-0091/ABC-0091_S6_L002_R2_001.fastq.gz
    local source_dir="$1"
    local destination_dir="$2"
    local file_type="$3"
    find "$source_dir" -name "${file_type}" | while read item; do
        parent_dir="${item##*Samples/}"
        parent_dir="${parent_dir%%/*}"
        output_path="${destination_dir}/${parent_dir}" 
        printf "Copying to:\n%s\n\n" "${output_path}"
        rsync  --dry-run -vctrPh "$item" "${output_path}/"  #
    done
}



# ~~~~~ GET SCRIPT ARGS ~~~~~ #
# project to be transfered from BaseSpace
project_ID="$1" # "NextSeq-123"

# directory which has been mounted by BaseMount (user should have done this already)
basepace_dir="$2"



# ~~~~~ MISC PARAMS ~~~~~ #
# type of file we want to copy
file_type='*.fastq.gz'



# ~~~~~ SETUP DIRS ~~~~~ #
# parent directory for sequencer output
sequencer_dir="/ifs/data/NGS_sequencer_dir"
# location where the project data will be copied to; make the dir
project_sequencer_dir="${sequencer_dir}/${project_ID}"
mkdir -p "$project_sequencer_dir"

# location where the BaseSpace projects will appear
basepace_projects="${basepace_dir}/Projects"

# BaseSpace dir for the project we want
project_BaseSpace_dir="${basepace_projects}/${project_ID}"



# ~~~~~ VALIDATIONS ~~~~~ #
# make sure project_BaseSpace_dir exists and contains projects
[ ! -d "$project_BaseSpace_dir" ] && printf "ERROR: Item is not a valid directory:\n%s\n\nDoes it exist? Exiting...\n\n" "$project_BaseSpace_dir" && exit
num_entries="$(ls -1 "$project_BaseSpace_dir" | wc -l)"
(( "$num_entries" < 1 )) && printf "ERROR: Directory contains %s items:\n%s\n\nExiting...\n\n" "$num_entries" "$project_BaseSpace_dir" && exit

# make sure that fastq files are present in the dir
num_fastqs="$(find "$project_BaseSpace_dir" -type f -name "*.fastq.gz" | wc -l)"
(( "$num_fastqs" < 1 )) && printf "ERROR: Directory contains %s fastq.gz files:\n%s\n\nExiting...\n\n" "$num_fastqs" "$project_BaseSpace_dir" && exit

# check to see if any filenames DO NOT contain the parent ID
# psuedo-boolean flag; "True" or "False"
some_filenames_lack_parentID="False"
echo "Checking to see if any fastq.gz files lack sample ID from their parent directory..."
find "$project_BaseSpace_dir" -type f -name "*.fastq.gz" | while read item; do
    # item="/ifs/home/kellys04/projects/Clinical_580_gene_panel/BaseSpace/Projects/NextSeq-123/Samples/ABC-0247/Files/ABC-0247_S12_L004_R2_001.fastq.gz"
    parent="${item##*Samples/}"
    parent="${parent%%/*}"
    file_base="$(basename "$item")"
    # echo "$parent"
    # echo "$file_base"
    # echo ""
    if echo "$file_base" | grep -q "${parent}"; then
        # echo "string ${parent} is in string $file_base"
        [ "$some_filenames_lack_parentID" != "True" ] && some_filenames_lack_parentID="True"
    fi
done

# the fastq files are currently in subdirs but we want to copy all of them to a common dir;
# check to make sure that none have duplicate names
echo "Checking to make sure that all fastq.gz files have unique names..."
num_all_fastq_files="$(find "${project_BaseSpace_dir}" -name "*.fastq.gz" | wc -l)"
num_unique_filenames="$(find "${project_BaseSpace_dir}" -name "*.fastq.gz" -exec basename {} \; | sort -u | wc -l)"



# ~~~~~ COPY THE FILES ~~~~~ #
# if False, then check for duplicate file basenames; if True, keep parent dirs
if [ "$some_filenames_lack_parentID" == "False" ]; then
    # if there are no duplicates, the values should be equal
    if [ "$num_all_fastq_files" -eq "$num_unique_filenames"  ]; then
        # it is safe to copy all fastq files to a common dir
        printf "Copying all %s files to common dir\n\n" "$file_type"
        copy_to_common_dir "$project_BaseSpace_dir" "$project_sequencer_dir" "$file_type"
    elif [ ! "$num_all_fastq_files" -eq "$num_unique_filenames"  ]; then
        # not equal; some fastq's have same filename, need to preserve parent dir's
        printf "Copying all %s files, preserving parent dirs\n\n" "$file_type"
        copy_with_parents "$project_BaseSpace_dir" "$project_sequencer_dir" "$file_type"
    fi
elif [ "$some_filenames_lack_parentID" == "True" ]; then
    printf "Copying all %s files, preserving parent dirs\n\n" "$file_type"
    copy_with_parents "$project_BaseSpace_dir" "$project_sequencer_dir" "$file_type"
fi

