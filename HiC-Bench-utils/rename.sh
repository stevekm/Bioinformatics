#!/bin/bash

# rename the dirs in the inputs dir

change_name () {
    local filename="$1"
    local old_pattern="$2"
    local new_pattern="$3"
    new_basename="$(echo "$(basename "$filename")" | sed -e "s|${old_pattern}|${new_pattern}|g")"
    new_filename="$(dirname "$filename")/${new_basename}"
    mv "$filename" "$new_filename"
}

rename_dir () {
    local old_pattern="$1"
    local new_pattern="$2"
    find fastq/ -maxdepth 1 -mindepth 1 -type d -name "*-${old}*" | while read item; do
        change_name "$item" "$old_pattern" "$new_pattern"
    done
}

old="K27AC" # fastq/PAUUDV-D-K9AC
new="H3K27AC" # fastq/PAUUDV-D-H3K27AC
rename_dir "$old" "$new"

old="K27ME3"
new="H3K27ME3"
rename_dir "$old" "$new"

old="K9AC"
new="H3K9AC"
rename_dir "$old" "$new"

old="CONT"
new="INPUT"
rename_dir "$old" "$new"

old="K27ME"
new="H3K27ME3"
rename_dir "$old" "$new"
