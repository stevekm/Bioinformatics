# Quality Control Scripts

Scripts for running Quality Control on NGS data in fastq.gz format

Current tools:
- FastQC: quality metrics for fastq
- Kraken: metagenomic analysis, useful for detecting possible contaminants
- Krona: create interactive pie chart from Kraken output

This script is designed to work on the NYU HPC system by submitting jobs to the SGE cluster, but it can also be modified to simply run in the current session. 

In the current implementation, the script will search the `input_dir` for all `*.fastq.gz` files, and run the selected script on them. A subdirectory for each input file will be create in the `output_dir`, where program output will be placed. 

Usage:
```bash
run_qc.sh /path/to/output_dir /path/to/input_dir FastQC_Kraken_Krona.sh
```

# Sample output

FastQC ouptput

![screen shot 2016-11-18 at 3 00 38 pm](https://cloud.githubusercontent.com/assets/10505524/20444653/e654f76a-ad9f-11e6-9e61-c49581d1a35c.png)

Krona plot snapshot

![rdp krona](https://cloud.githubusercontent.com/assets/10505524/20355138/57fb07be-abee-11e6-8c88-18ba2344d7d0.png)
