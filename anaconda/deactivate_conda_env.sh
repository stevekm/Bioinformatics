#/bin/bash

# `source` this to deactivate the conda environment you started with the other script

source deactivate "$CONDA_ENV_NAME"

export PATH="$OLD_PATH"
export PYTHONPATH="$OLD_PYTHONPATH"

python_module="python/2.7"
module load "$python_module"
