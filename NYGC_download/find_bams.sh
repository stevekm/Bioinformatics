#!/bin/bash

(
ssh username@scp.nygenome.org <<E0F
set -x
find /data/delivery/NYGC/data/ -name "*.bam" -type f -exec readlink -f {} \;
E0F
) > find_bam.txt
