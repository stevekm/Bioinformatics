Benchmarking script for `GetBaseCountsMultiSample` program

https://github.com/zengzheng123/GetBaseCountsMultiSample

Put a list of bam, bai, and sample ID's into `all_bam_samples.tsv`, and put your variants to benchmark with in `consensus.maf`

Setup with `make install singularity-pull`, run with `make run`, and plot with `./plot_trace.R trace.txt`


- https://hub.docker.com/r/cmopipeline/getbasecountsmultisample/tags

- https://hub.docker.com/layers/mskcc/helix_filters_01/getbasecountsmultisample-1.2.2/images/sha256-32c64ff26ce85c791cadb597a2b665ef516e396ace710ffe3cdc155b95c9d014?context=repo


## Results

Tested with maf file containing 417632 variants, and 48 targeted exome sequencing .bam files ranging from 600MB - 8GB in size

- `aggregate.tsv`: metrics showing mean execution time per sample

- `total_time.pdf`: total execution time for HPC jobs with varying numbers of samples and CPU threads

- `time_per_sample.pdf`: mean execution time per sample
