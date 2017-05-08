#!/bin/bash

# this script will find the peaks files in the given directory and run the annotator on them

# ~~~~~~~~~~ CUSTOM FUNCTIONS ~~~~~~~~~~ #
error_on_zerolength () {
    local test_string="$1"
    [ -z "$test_string" ] && echo -e "ERROR: String is length zero\nExiting..." && exit
}


file_exists () {
    local input_file="$1"
    error_on_zerolength "$input_file"
    [ ! -f "$input_file" ] && echo -e "ERROR: Item is not a file:\n$input_file\nDoes it exist?\nExiting..." && exit
}

dir_exists () {
    local input_dir="$1"
    error_on_zerolength "$input_dir"
    [ ! -d "$input_dir" ] && echo -e "ERROR: Item is not a dir:\n$input_dir\nDoes it exist?\nExiting..." && exit
}


annotate_beds () {
    local input_dir="$1"
    find "$input_dir" -name "*.bed" -print0 | while read -d $'\0' item; do
        item_basename="$(basename "$item")"
        output_name="$(dirname "$item")/${item_basename%%.bed}_annotation.tsv"
        $annotate_script "$item" "$output_name"
    done
}


# ~~~~~ CHECK SCRIPT ARGS ~~~~~ #
num_args="$#"
args_should_be_greaterthan="0"
if (( "$num_args" <= "$args_should_be_greaterthan" )); then
            echo "ERROR: Wrong number of arguments supplied"
            echo "Number of script arguments should be at least: $args_should_be_greaterthan"
            grep '^##' $0
            exit
fi

# ~~~~~ GET SCRIPT ARGS & SETTINGS ~~~~~ #
input_dirs="${@:1}" # accept a space separated list of dirs
annotate_script="$(readlink -f annotate_peaks.R)"
file_exists "$annotate_script"

# ~~~~~~~~~~ RUN ~~~~~~~~~~ #
for input_dir in $input_dirs; do
    dir_exists "$input_dir"
    annotate_beds "$input_dir"
done
