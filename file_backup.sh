#!/bin/bash


## USAGE: file_backup.sh <file_to_backup> <backup_output_dirname>


# ~~~~~~ script args processing ~~~~~~ #
# # if not enough args, output USAGE info lines (which start with '##') and exit
if (($# != 2)); then
  grep '^##' $0
  exit
fi


file_to_backup="$1"
backup_dir="$2"

mkdir -p "$backup_dir"

# get the basename of the file to backup
basename_old_file="$(basename $file_to_backup)"
#echo "$basename_old_file"

# append the timestamp
basename_new_file="${basename_old_file}_$(date -u +%Y%m%dt%H%M%S)"
#echo "$basename_new_file"

# move the old file to the backup dir
echo -e "Moving $file_to_backup to ${backup_dir}/${basename_new_file}"
mv "$file_to_backup" "${backup_dir}/${basename_new_file}" && echo "File backed up to: ${backup_dir}/${basename_new_file}"
