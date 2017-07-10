#!/bin/bash
set -x

# the command to use for transfer

# first get everything except the huge BAM files
# rsync -zrtvhP -e ssh username@scp.nygenome.org:/data/delivery/SomeFoundation/data/Project_123 . --exclude="*.bam"

# then get the rest of the files
rsync -zrtvhP -e ssh username@scp.nygenome.org:/data/delivery/SomeFoundation/data/Project_123 .
