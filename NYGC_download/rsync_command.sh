#!/bin/bash
set -x

# the command to use for transfer
rsync -zrtvhP -e ssh username@scp.nygenome.org:/data/delivery/SomeFoundation/data/Project_123 .
