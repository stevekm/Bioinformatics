---
title: "Transfer Files from Basespace"
author: "Stephen Kelly"
date: "3/28/2017"
output: html_document
---

This workflow will show you how to mount your BaseSpace account on the HPC server so you can copy files (Projects and Runs) from your BaseSpace account to the server storage for processing. 

BaseMount Help Page:

https://help.basespace.illumina.com/articles/descriptive/introduction-to-basemount/

Before you start:
- log into your Basespace account in a web browser: https://basespace.illumina.com/home/index
- make sure your account has loaded the associated group that has your desired projects
- - e.g. if your project is under your Personal account, then switch to your Personal account in BaseSpace before you start

# Loading BaseMount

BaseMount is already installed on the phoenix server. You just need to create and mount a directory to load the files into 

```bash
project_dir="/ifs/home/kellys04/projects/sequencing"
basepace_dir="${project_dir}/BaseSpace"
mkdir -p "$basepace_dir"
basemount "${basepace_dir}/"
```

Now you will get a BaseMount message with a URL to open in your BaseSpace account; copy/paste the URL into the web browser

Should look like this:

```bash

,-----.                        ,--.   ,--.                         ,--.
|  |) /_  ,--,--. ,---.  ,---. |   `.'   | ,---. ,--.,--.,--,--, ,-'  '-.
|  .-.  \' ,-.  |(  .-' | .-. :|  |'.'|  || .-. ||  ||  ||      \'-.  .-'
|  '--' /\ '-'  |.-'  `)\   --.|  |   |  |' '-' ''  ''  '|  ||  |  |  |
`------'  `--`--'`----'  `----'`--'   `--' `---'  `----' `--''--'  `--'
Illumina BaseMount v0.15.15.1872 public  2016-12-16 10:47

Command called:
    basemount /ifs/home/kellys04/projects/sequencing/BaseSpace/
From:
    /ifs/home/kellys04/projects/sequencing

Starting authentication.

You need to authenticate by opening this URL in a browser:
  https://basespace.illumina.com/oauth/device?code=xxxx
........
It worked!
Your identification has been saved.

Mounting BaseSpace account.
To unmount, run: basemount --unmount /ifs/home/kellys04/projects/sequencing/BaseSpace

```

# Copy the Files

Make a list of the fastq files to get for the Project you want, and copy them with `rsync` to the sequencer directory.

Since we want to preserve the directory structure from BaseSpace, using the Project ID as the parent directory, we will first assemble a list of the fastq files to be copied with relative paths. Then we will run `rsync` from within the mounted BaseSpace directory to copy the files. 

Its strongly advised to run this from within `screen` on the server. 

```bash
sequencer_dir="/ifs/data/mylab/NexSeq2"

basepace_dir="/ifs/home/kellys04/projects/sequencing/BaseSpace"
basepace_projects="/ifs/home/kellys04/projects/sequencing/BaseSpace/Projects"
files_list_dir="/ifs/home/kellys04/projects/sequencing/BaseSpace_files_lists"; mkdir -p "$files_list_dir"

project_ID="NGS-Project-1"
project_files_list="${files_list_dir}/${project_ID}.txt"


(
cd "${basepace_projects}"

# make list of fastq files for the run; preserve the 
find "${basepace_projects}/${project_ID}" -name "*.fastq.gz" | sed -e "s|${basepace_projects}/||g" > "$project_files_list"

printf "\nFiles to be copied:\n\n"
cat "$project_files_list"

# transfer the files
rsync --progress -vRt --files-from="$project_files_list" "$basepace_projects" "$sequencer_dir" # --dry-run
)


```

# Unmount BaseSpace

When you are done copying, you need to unmount the BaseSpace directory. Use this command

```
basemount --unmount /ifs/home/kellys04/projects/sequencing/BaseSpace
```

# Notes

Authentication

The first time you run BaseMount, you will be directed to a web URL and asked to enter your BaseSpace Sequence Hub user credentials. BaseMount will use these credentials to authenticate your interactions with BaseSpace Sequence Hub. By default, the credentials are cached in your home directory and they can be password-encrypted for security, just like an ssh key.

- this means that if you need to reset your login credentials on the server with BaseMount, you can just nuke the `~/.basespace` dir. It will get re-created when you start the procedure over again. 
