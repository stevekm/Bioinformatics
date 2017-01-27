#!/bin/bash

# manual: http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html

# here are some notes about using heredocs and qsub and Rscript

# you can submit qsub jobs directly from a heredoc:

my_var="foobar"
qsub <<HERE
  #!/bin/bash
  echo "$my_var"

HERE

# you can run R code from a heredoc within your bash script

my_var="foobar"
Rscript - "arg1" "arg2" <<E0F
  ## R code
  cat("\nR loaded\n")
  args <- commandArgs(TRUE)
  cat("Script args are:\n")
  args
  print("$my_var")
  sessionInfo()
E0F

# btw here are some other ways to pass heredoc to R
# R --slave <<EOF
# Rscript --slave --vanilla - <<EOF

# you can can submit an R script to run from qsub as well:

qsub -b y -wd "$tmp_outdir" -o :${tmp_logdir}/ -e :${tmp_logdir}/ -pe threaded 1 Rscript my_script.R "$arg_1" "$arg_2"

# NOTE: '-b y' is required here
# also your R script needs to start like this:
#!/usr/bin/env Rscript


# things I haven't tested yet but will probably need eventually:

# nested heredocs:
# http://stackoverflow.com/questions/13254077/is-nesting-of-here-document-possible-in-a-unix-bash-script
ssh s1 << \eof1
  ssh s2 << \eof2
    hostname
eof2
eof1

# and disable parameter expansion within the heredoc by quoting the initial ... thing
# this will probably be useful for R code where the '$' is used, since that breaks in heredocs
# http://tldp.org/LDP/abs/html/here-docs.html
cat <<'Endofmessage'

Hello, there, $NAME.
Greetings to you, $NAME, from $RESPONDENT.

Endofmessage

# TO-DO: try to qsub an R script heredoc directly.. 
