#/bin/bash

# `source` this to start the conda environment

env_name="${1:-python3.6}"

export CONDA_ENV_NAME="$env_name"
export OLD_PATH="$PATH"
export OLD_PYTHONPATH="$PYTHONPATH"

module unload python
export PYTHONPATH=/ifs/home/kellys04/anaconda3/bin:$PYTHONPATH
export PATH=/ifs/home/kellys04/anaconda3/bin:$PATH

# conda create -n myenv python=3.6
source activate "$env_name"
