#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script will run Delly2 for fusion detection
on a directory of .bam files
and convert the output .bcf files to .vcf with bcftools

Examples
--------
Example usage::

    $ ./run_Delly2.py BAM-STAR

Where ``BAM-STAR`` is a directory containing .bam files with their corresponding .bai files
"""
import sys
import os
import subprocess as sp
from time import sleep

# ~~~~ GLOBALS ~~~~~~ #
config = {
'delly2_bin': '/ifs/home/kellys04/software/delly/src/delly',
'bcftools_bin': '/ifs/home/kellys04/software/delly/src/bcftools/bcftools',
'hg19_fa': '/local/data/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa',
'call_types': (
    ('deletions', 'DEL'),
    ('duplications', 'DUP'),
    ('inversions', 'INV'),
    ('translocations', 'BND'),
    ('insertions', 'INS')
)
}



# ~~~~ FUNCTIONS ~~~~~~ #
def subprocess_cmd(command):
    """
    Runs a subprocess command

    Parameters
    ----------
    command: str
        a shell command to run

    Returns
    -------
    tuple
        a tuple containing the stdout and stderr messages in the format ``('stdout', 'stderr')``
    """
    process = sp.Popen(command, stdout = sp.PIPE, stderr = sp.PIPE, shell = True, universal_newlines = True)
    proc_stdout, proc_stderr = process.communicate()
    return(proc_stdout, proc_stderr)

def submit(command, log_dir, job_name = 'python', params = '-wd $PWD -j y', pre_commands = 'set -x', post_commands = ''):
    """
    Submits a qsub job to run on the HPC cluster (SGE)

    Parameters
    ----------
    command: str
        the shell command to be run in the qsub job
    log_dir: str
        path to directory to hold qsub logs
    job_name: str
        the name of the qsub job; must not start with a number!
    params: str
        extra params to use with ``qsub``
    pre_commands: str
        commands to prepend to the start of the job commands
    post_commands: str
        commands to append to the end of the job commands
    """
    stdout_log_dir = log_dir
    stderr_log_dir = log_dir


    qsub_command = """
qsub {0} -N "{1}" -o :"{2}" -e :"{3}" <<E0F
{4}
{5}
{6}
E0F
""".format(
params,  # 0
job_name, # 1
stdout_log_dir, # 2
stderr_log_dir, # 3
pre_commands, # 4
command, # 5
post_commands # 6
)

    print(qsub_command + '\n')
    proc_stdout, proc_stderr = subprocess_cmd(command = qsub_command)
    print((proc_stdout, proc_stderr))
    print('\n\n')
    sleep(1)

def mkdirs(path, return_path=False):
    """
    Make a directory, and all parent dir's in the path

    Parameters
    ----------
    path: str
        the path to a directory to create
    return_path: bool
        whether or not the created path should be ``return``'ed

    Returns
    -------
    str
        the path that was created, if ``return_path`` was passed as ``True``
    """
    import sys
    import os
    import errno
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
    if return_path:
        return path

def find_bams(bam_dir):
    """
    Finds all the .bam files in the directory

    Parameters
    ----------
    bam_dir: str
        path to directory containing .bam files. Should also contain corresponding .bai files

    Returns
    -------
    list
        a list of the paths to all the .bam files
    """
    bam_files = []
    for root, dirs, files in os.walk(bam_dir):
        for file in files:
            if file.endswith(".bam"):
                bam_files.append(os.path.join(root, file))
    return(bam_files)

def delly2_cmd(bam_file, output_dir, call_type_name, call_type_arg, sampleID = None):
    """
    Build the terminal commands to run Delly2 on a single sample
    Make a separate command for each SV calling type, concatenate them all together

    Parameters
    ----------
    bam_file: str
        path to a .bam file
    output_dir: str
        path to an output_dir
    call_type_name: str
        the name of the call typ to use, e.g. 'deletions'
    call_type_arg: str
        the CLI arg to use for the command, e.g. 'DEL'
    sampleID: str
        identifier for the sample

    Returns
    -------
    str
        a shell command to run Delly2, followed by bcftools, for the .bam file
    """
    # get params from config
    delly2_bin = config['delly2_bin']
    bcftools_bin = config['bcftools_bin']
    hg19_fa = config['hg19_fa']

    output_SV_bcf_ext = '.bcf'
    output_SV_vcf_ext = '.vcf'

    # get the sampleID from the bam_file if not passed
    if not sampleID:
        sampleID = os.path.splitext(os.path.basename(bam_file))[0]

    # make a command for each SV calling type:

    # make filepath for the .bcf output
    sample_output_SV_bcf_basename = ''.join([sampleID, '.' + call_type_name, output_SV_bcf_ext])
    sample_output_SV_bcf = os.path.join(output_dir, sample_output_SV_bcf_basename)

    # make filepath for the .vcf converted output
    sample_output_SV_vcf_basename = ''.join([sampleID, '.' + call_type_name, output_SV_vcf_ext])
    sample_output_SV_vcf = os.path.join(output_dir, sample_output_SV_vcf_basename)

    # delly call -t DEL -g "genome.fa" -o "results_dir/delly2-snv/Sample1.deletions.bcf" "results_dir/BAM-GATK-RA-RC/Sample1.dd.ra.rc.bam"
    # bcftools view "results_dir/delly2-snv/Sample1.deletions.bcf" > "results_dir/delly2-snv/Sample1.deletions.vcf"
    command = """
{0} call -t {1} -g "{2}" -o "{3}" "{4}"
{5} view "{6}" "{7}"
""".format(
delly2_bin,
call_type_arg,
hg19_fa,
sample_output_SV_bcf,
bam_file,

bcftools_bin,
sample_output_SV_bcf,
sample_output_SV_vcf
)
    return(command)

def make_Delly2_cmds(bam_files, output_dir):
    """
    Makes the Delly2 commands for all the bam files

    Parameters
    ----------
    bam_files: list
        a list of paths to .bam files
    output_dir: str
        a path to an output directory

    Returns
    -------
    dict
        a dictionary with entries in the format of ``'job_name': 'command'``

    """
    call_types = config['call_types']
    commands = {}
    for bam_file in bam_files:
        for call_type_name, call_type_arg in call_types:
            sampleID = os.path.splitext(os.path.basename(bam_file))[0]
            job_name = 'Delly2_{0}_{1}'.format(sampleID, call_type_name)
            commands[job_name] = delly2_cmd(bam_file = bam_file,
                                            output_dir = output_dir,
                                            call_type_name = call_type_name,
                                            call_type_arg = call_type_arg,
                                            sampleID = sampleID)
    return(commands)




def main(**kwargs):
    """
    Main control function for the program. Submits a qsub job to the HPC to run Delly2 & bcftools on each sample

    Keyword Args
    ------------
    bam_dir: str
        path to directory containing .bam files
    output_dir: str
        path to directory to hold output
    log_dir: str
        path to directory to hold qsub logs


    Notes
    -----
    Creates and runs a shell command that looks like this::

        qsub -wd $PWD -j y -N "Delly2_174_inversions" -o :"/ifs/data/molecpathlab/PNET_GYN/Delly2_fusions/logs" -e :"/ifs/data/molecpathlab/PNET_GYN/Delly2_fusions/logs" <<E0F
        set -x

        /ifs/home/kellys04/software/delly/src/delly call -t INV -g "/local/data/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa" -o "/ifs/data/molecpathlab/PNET_GYN/Delly2_fusions/output/174.inversions.bcf" "BAM-STAR/174.bam"
        /ifs/home/kellys04/software/delly/src/bcftools/bcftools view "/ifs/data/molecpathlab/PNET_GYN/Delly2_fusions/output/174.inversions.bcf" "/ifs/data/molecpathlab/PNET_GYN/Delly2_fusions/output/174.inversions.vcf"


        E0F


        ('Your job 4133697 ("Delly2_174_inversions") has been submitted\n', '')

    """
    # get the args that were passed
    bam_dir = kwargs.pop('bam_dir', None)
    output_dir = kwargs.pop('output_dir', "output")
    log_dir = kwargs.pop('log_dir', "logs")

    if not bam_dir:
        print('no bam_dir provided')
        sys.exit()
    if not os.path.isdir(bam_dir):
        print('bam_dir is not a dir: {0}'.format(bam_dir))
        sys.exit()

    # make the other dir paths if they dont exist
    for path in [output_dir, log_dir]:
        mkdirs(path)

    # reset for the full paths
    output_dir = os.path.realpath(os.path.expanduser(output_dir))
    log_dir = os.path.realpath(os.path.expanduser(log_dir))

    # get the .bam files from the dir
    bam_files = find_bams(bam_dir = bam_dir)

    # make the Delly2 shell commands for each file
    Delly2_cmds = make_Delly2_cmds(bam_files = bam_files, output_dir = output_dir)

    # run all the jobs
    for job_name, command in Delly2_cmds.items():
        submit(command = command,
                log_dir = log_dir,
                job_name = job_name)





def parse():
    """
    Parses CLI arguments and passes them to ``main()``
    """
    bam_dir = sys.argv[1]
    output_dir = "output"
    log_dir = "logs"

    kwargs = {
    'bam_dir': bam_dir,
    'output_dir': output_dir,
    'log_dir': log_dir
    }
    main(**kwargs)

# ~~~~ RUN ~~~~~~ #
if __name__ == "__main__":
    parse()
