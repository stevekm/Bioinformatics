#!/usr/bin/env python
# python 2.7

import sys
import os
import argparse


# read the lines from the heredoc passed on stdin
stdin_items = []
for line in sys.stdin:
    stdin_items.append(line)

# remove the '-' arg from sys.argv since it represents the heredoc we passed
sys.argv =[x for x in sys.argv if x != '-']

# parse the rest of the args however you want
parser = argparse.ArgumentParser(description='Demo script to read stuff from heredoc along with other args')
parser.add_argument("file_list", help="Paths to input files", nargs="+")
args = parser.parse_args()
file_list = args.file_list


print "File list items:"
for item in file_list:
    print item


print "Items read from stdin heredoc:"
for item in stdin_items:
    print item


