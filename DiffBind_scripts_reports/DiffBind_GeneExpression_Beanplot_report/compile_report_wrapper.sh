#!/bin/bash
# set -x

##
## USAGE: compile_report_wrapper.sh document.Rnw
## This script will compile a .Rnw document with knitr in R, then compile the results .tex file with pdflatex
##

# check script args
if [ $# != 1 ] # if not enough args provided
then
  echo -e "Only one argument is supported for this script at this time" # update this in this in the future, use like '$@' in a for loop
  grep '^##' $0 # print out lines from the script that start with '##' and exit
  exit
fi

# get script arg
tmp_file="$1"


# compile the document with knitr
# in the future update to make sure that the arg passed was a file that ends in .Rnw
Rscript --slave --no-save --no-restore - "$1" <<EOF
  ## R code
  cat("\nR loaded\n")
  args <- commandArgs(TRUE)
  cat("Script args are:\n")
  args
  cat("Checking to see if R pacakge knitr is available...\n")
  # check if knitr is installed.. 
  # if (!require("knitr")) install.packages("knitr")
  if (!require("knitr")) {
    # cat("knitr not installed, attempting to install..\n\n")
    cat("Sorry you've got to go install knitr to your system before this script can compile the document\n\n")
    # code that I tried to use to install R locally.. this doesn't work... 
    # system("mkdir -p ~/software/R/R_libs")
    # tmp_dir<-system("readlink -m ~/software/R/R_libs",intern = TRUE)
    # cat("\n\n")
    
    # install.packages("knitr",lib = tmp_dir,repos="http://cran.rstudio.com/",destdir = tmp_dir)
    # library("knitr",lib.loc = tmp_dir)
    
    # exit if knitr not installed, haven't fogured out yet how to get the script to install it locally for the user / session
    quit()
  } else {
  cat("knitr is already installed, loading knitr\n")
  library("knitr")
    
  }
  cat("Compiling document with knitr\n")
  knit(args[1])
EOF


# btw here are some other ways to pass heredoc to R
# R --slave <<EOF
# Rscript --slave --vanilla - <<EOF


# PDF compilation with pdflatex
# check if the .tex file exists; 
# # sometimes knitr can fail but still produce a tex file
if [ -f ${tmp_file%%.Rnw}.tex ]; then
  # check if pdflatex is installed
  if command -v pdflatex &>/dev/null; then
      # gdate "$@"
      echo -e "Compiling the TeX file; file is:\n${tmp_file%%.Rnw}.tex\n\n"
      # compile with pdflatex
      pdflatex "${tmp_file%%.Rnw}.tex"
      # compile it again...
      pdflatex "${tmp_file%%.Rnw}.tex"
      # one more time just to be safe... yes, this matters!
      pdflatex "${tmp_file%%.Rnw}.tex"
  else
      # date "$@"
      echo -e "pdflatex not installed, unable to compile file!\n"
  fi
else
  echo -e "TeX file not found! File should be:\n${tmp_file%%.Rnw}.tex\n\n"
fi
