#!/usr/bin/env Rscript

## USAGE: compile_RMD_report.R report.Rmd
# module load pandoc/1.13.1

# ~~~~~ GET SCRIPT ARGS ~~~~~~~ #
args <- commandArgs(TRUE)

Rmdfile <- args[1]

rmarkdown::render(Rmdfile)