## HOMER_motif_analysis
- A collection of scripts to analyze peaks from a .bed file; run HOMER `findMotifsGenome.pl` to discover motifs, then aggregate the discovered motifs and annotate their genomic locations with `annotatePeaks.pl`

- These scripts are designed to run on a HPC cluster using a variant of SGE job scheduler with `qsub`

- Set up everything in the `workflow_scratchpad.sh` file, then run all the commands on the shell. This is nice because it keeps a verbose record of the relevant shell settings prior to running the scripts, including the `qsub` commands and settings. 

- `motif_analysis_HOMER2_params.sh` contains HOMER settings and a few other things
 
- `motif_analysis_HOMER3.sh` is the script that does all the analysis and gets submitted to `qsub` for completion on the HPC
