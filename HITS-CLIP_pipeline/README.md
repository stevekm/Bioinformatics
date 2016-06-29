##HITS-CLIP pipeline

####Files:

- `workflow.Rmd` : A simple text file where I hard-code in the locations of the directories and scripts needed, and then set up needed code to loop over all samples and create create output directories and submit `qsub` jobs to the cluster. 

- `samplesheet.tsv` : (not included) A tab-delimited file containing one row per sample, with relevant information such as location to reference databases specific to the sample, genome version, fastq file locations, etc., which gets parsed out and passed to the pipeline script.

- `HITSCLIP_pipeline_trim16_qsub.sh` : A bash script to be submitted to the cluster, which recieves arguments parsed out by the `workflow` scripts and executes the pipeline. 



####Usage:

- Set up and modify `samplesheet` and `workflow` as needed to group the samples and create the directory structure.

- Run the `workflow` to submit `pipeline` script jobs for running on the cluster.



####Pipeline Overview:

(modify these steps as necessary)

- run FastQC on the input fastq files

- filter out low quality reads, collapse duplicates, and trim reads as necesary

- run Kraken metagenomics analysis on the fastq files to analyze for potential contaminants

- (more stuff, I'll fill this in later.. standard HITS-CLIP pipeline from publication)

- Run MEME-ChIP & HOMER motif analysis


