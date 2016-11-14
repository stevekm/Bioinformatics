#!/bin/bash
set -x 
## USAGE: run.sh /path/to/output_dir /path/to/input_dir /path/to/runscript.sh
# this script might be broken, fix it later.. 
# this script will loop over all fastq.gz files it finds and run the script
# putting the results in an outdir with the same name as the input file

OUTDIR="$1"
mkdir -p "$OUTDIR"

INPUTDIR="$2"
FILES="$(find $INPUTDIR -name "*.fastq.gz")"

# tmp_script="$(dirname $0)/FastQC_Kraken_qsub.sh"
tmp_script="$3"
chmod +x "$tmp_script"
cat $tmp_script

for i in $FILES; do
    tmp_fastq="$(readlink -f $i)"
    tmp_outdir="${OUTDIR}/$(basename $tmp_fastq)"
    tmp_logdir="${tmp_outdir}/logs"
    mkdir -p "$tmp_logdir"

    qsub -wd $tmp_outdir -o :${tmp_logdir}/ -e :${tmp_logdir}/ -pe threaded 2-4 "$tmp_script" "$tmp_outdir" "$tmp_fastq"
    
    # run without qsub here:
    # "$tmp_script" "$tmp_outdir" "$tmp_fastq" 
done
