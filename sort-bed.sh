#!/bin/bash

sort_bed () {
    local bed_file="$1"
    sort -k1,1 -k2,2n "$bed_file" > tmp && mv tmp "$bed_file"
}

bed_file="$1"
sort_bed "$bed_file"
