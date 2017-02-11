---
title: "Workflow"
author: "Stephen Kelly"
date: "8/29/2016"
output: html_document
---

Dual barplots for alignment sample stats vs. peaks per sample. Code for multiple types of barplots included. 

Usage:
```
code/peaks_alignment_barplots_workflow.R input/peaks_stats.tsv input/alignment_stats_extended.csv  SamplePeaks output/

```

Requires formatted input files (see examples) for alignment stats and peak stats


To use this with HiC-Bench:

- copy this directory
- symlink to the parent align-stats results dir
```bash
align-stats-results -> ../../pipeline/align-stats/results_OG_good/align-stats.standard/align.by_sample.bowtie2
```
- symlink to the parent peaks results dir
```bash
peaks-results -> ../../pipeline/peaks/results/peaks.by_sample.macs_broad/align.by_sample.bowtie2
```
- make peaks stats file
```bash
$ echo -e 'Peaks\tSample' > input/peaks_stats.tsv
$ find peaks-results/ -type f -name "peaks.bed" -exec wc -l {} \; | sed -e 's|peaks-results/||g' -e 's|/peaks.bed||g' | tr ' ' '\t' >> input/peaks_stats.tsv
```
