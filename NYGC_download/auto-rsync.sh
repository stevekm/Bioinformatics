#!/bin/bash
set -x

RSYNC="rsync -zrtvhP -e ssh username@scp.nygenome.org:/data/delivery/SomeFoundation/data/Project_123 ."

echo "Starting first attempt"
$RSYNC
RV=$? # exit status of first attempt
CNT=1
if [ $RV -eq 0 ]; then
    echo "Finished successfully on first pass."
else
    while [ $RV -ne 0 ]; do
        CNT=$(( CNT + 1 ))
        echo "Starting attempt $CNT"
        $RSYNC
        RV=$?
    done
    echo "Finished w/ zero return value on attempt $CNT"
fi
