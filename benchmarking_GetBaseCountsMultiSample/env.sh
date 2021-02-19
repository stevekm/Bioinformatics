#!/bin/bash
# environment needed to run benchmarking on Juno / Silo server
if module avail java/jdk1.8.0_202 1&>/dev/null; then module load java/jdk1.8.0_202; fi
module load singularity/3.3.0
