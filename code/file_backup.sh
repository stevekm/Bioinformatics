#!/bin/bash


## USAGE: file_backup.sh <file_to_backup> <backup_output_dirname>
## EXAMPLE: file_backup.sh file_to_keep.txt old
## DESCRIPTION: this is Steve's handy script for stashing old files that you don't want to delete because he's a packrat that hates deleting files they might be useful again someday
## backup dir is automatically created and file_to_keep is moved there with a timestamped filename

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

# get the extension of the old file
old_ext="${basename_old_file##*.}"

# append the timestamp
basename_new_file="${basename_old_file}_$(date -u +%Y%m%dt%H%M%S).${old_ext}"
#echo "$basename_new_file"

# move the old file to the backup dir
printf "Moving $file_to_backup to %s/%s" "${backup_dir}" "${basename_new_file}"
mv "$file_to_backup" "${backup_dir}/${basename_new_file}" && echo "File backed up to: ${backup_dir}/${basename_new_file}"
