Examples of workflows to run tools on the cluster, with implementation notes. 

## Motif Analysis with HOMER

```bash
#!/bin/bash

# Motif analysis workflow for Super Enhancer BED file inputs 

# directory we're doing the work in # (1.)
project_dir="/ifs/home/kellys04/projects/SmithLab_ChIpSeq_2016-12-31/project_notes/motifs-SuperEnhancers"

# directory that holds our input BED files # (1.)
SE_dir="/results/ChIP-Seq/superenhancer/ROSE/bed/"

output_HOMER="${project_dir}/HOMER"
mkdir -p "$output_HOMER" # (2.)

# (3.)
run_HOMER_motif () {
    # function to run HOMER on a BED file, and set up output directories
    # http://homer.salk.edu/homer/ngs/peakMotifs.html
    local input_bed="$1" # (4.)
    local genome="hg19"
    local motif_options="-size given"
    local sample_ID="$(basename "$input_bed" | cut -d '-' -f1)"
    local sample_status="$(basename "$input_bed" | cut -d '-' -f2)"
    local outdir="${output_HOMER}/${sample_ID}/${sample_status}"
    mkdir -p "$outdir"
    local preparsed_dir="${outdir}/preparsed"
    mkdir -p "$preparsed_dir"
    local logdir="${outdir}/logs"
    mkdir -p "$logdir"
    local job_name="homer-${sample_ID}"
    echo "$sample_ID $sample_status"
    echo "$outdir"
    # (5.)
    qsub -wd $outdir -o :${logdir}/ -e :${logdir}/ -j y -pe threaded 4 -N "$job_name" <<E0F
# (6.)
module load homer/v4.6 
set -x
threads=\$NSLOTS 
findMotifsGenome.pl $input_bed $genome $outdir/ -p $threads -preparse -preparsedDir ${preparsed_dir}/ -dumpFasta $motif_options # && rm -rf $preparsed_dir
set +x
E0F
}


# we need to find all the input BED files, and pass them to our function to submit a qsub job for HOMER on it
# file paths look like this:
# /results/ChIP-Seq/superenhancer/ROSE/bed/ABC-R-H3K27AC_ROSE_superenhancer.bed

# HOMER
# (7.)
find "$SE_dir" -name "*.bed" -exec readlink -f {} \; | while read item; do
    echo "$item" 
    run_HOMER_motif "$item"
done
```
Notes:
1. Its a good idea to start with the full path to the directories you're working on. 
2. Use `mkdir -p` to create all directories in the given path.
3. Isolate the code needed to perform your task within a function, to make the code more readable and keep 'file-finding' code from getting tangled up in your 'task-performing' code. Also note that we're using the "old style" [bash function syntax](http://tldp.org/LDP/abs/html/functions.html) here, since it is more portable in case you end up working on an older system. 
4. Inside the function, use the `local` command to set variables in a scope local to the function, so they are not globally accessible by the rest of your code. This way, if something goes wrong while your code is running, you don't accidentally use `input_bed` that was set by another iteration of the function.
5. Command to submit the job to `qsub` to run on the compute cluster. 
    1. `-wd`: Sets the working directory for the job to run in. By default SGE runs jobs out of a type of temp directory, which can cause issues with some programs. You can set the working directory for the job explicitly with an argument like `-wd $outdir`, or you can use `-cwd` to run from the current working directory
    2. `-o :${logdir}/` sets the location for the std out log (and `-e` for std error). The trailing slash and leading colon are required. 
    3. `-j y`: merges the std out and std error logs into a single file (so `-e :${logdir}/` was not actually required here)
    4. `-pe threaded 4`: The number of threads for the job to use, can also be a range such as `4-8`
    5. `-N "$job_name"`: The name of the job to use, instead of the default value. This is the name that appears in `qstat`
6. We are using a [bash "heredoc"](http://tldp.org/LDP/abs/html/here-docs.html) to pass the commands for `qsub` instead of using an external script file. This is useful since it keeps all your code in one place, and allows you to dynamically build your commands to be run. Heredoc's do have a few quirks to be aware of if you plan to use them, so be sure to look over the [documentation](http://tldp.org/LDP/abs/html/here-docs.html) for them.
    1. The code will be running in a new bash session. First we need to load the `module` for our program, and set other environment variables
    2. Using `set -x` is useful for debugging since it prints all commands run by bash to std out, with all variables filled in. These entries are prefixed by `+` in the console or log file. This way, when you have a lot of dynamically set variables, you can see exactly what is getting run. You turn this feature back off with `set +x`. 
    3. By default, the heredoc will try to fill in any items prefixed with `$` as if it were a variable set in the session in which it was created. So, `$input_bed` will get filled in with the value that we assigned previously. This feature can be disabled on a per-item basis by escaping with a `\` like this: `\$`
    2. `$NSLOTS` is an environment variable set by `qsub` within the running job, and represents the number of threads ("slots") assigned to the job; we want to assign this value to the `threads` variable after the heredoc has been submitted as a job, so it is escaped and written as `threads=\$NSLOTS`. If you wanted to have a default value for `threads` to use with the code when `NSLOTS` isn't set or doesn't exist, you could use `threads=${NSLOTS:=4}`, which in this case would also need to be escaped for our heredoc. 
    3. As described, all the variables passed to the `findMotifsGenome.pl` program here will be filled in either from the function, or from variables assigned within the heredoc.
7. Here we use a simple `find` command to search for all the BED files in the specified directory. 
    1. The argument `-exec readlink -f {} \;` is used to ensure that the full path is returned for every file found. Not necessarily required in this example since our `SE_dir` variable started with the full path in the first place, but this is a good practice to use when passing around file paths to be processed.
    2. Every line output by `find` is passed to a `while` loop, which will print the file being processed and then run the `run_HOMER_motif` function on the file. Using `-exec run_HOMER_motif` is not possible here since `-exec` only works with programs, commands, and script invocations that are accessible external to the current session; we're trying to keep all our code in one place, so we want to avoid that.  
