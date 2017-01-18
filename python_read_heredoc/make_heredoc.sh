#!/bin/bash

program1_args="-d thing.txt -s foobar"
program2_args="-s another_arg -i some_ID"

py_script="$(dirname $0)/read_heredoc.py"

$py_script foo.txt bar.txt - <<E0F
$program1_args
$program2_args
E0F
