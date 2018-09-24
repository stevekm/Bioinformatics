SHELL:=/bin/bash
NXF_VER:=0.31.1
UNAME:=$(shell uname)
none:


# ~~~~~ SETUP CONDA ~~~~~ #
ifeq ($(UNAME), Darwin)
CONDASH:=Miniconda2-4.5.4-MacOSX-x86_64.sh
endif

ifeq ($(UNAME), Linux)
CONDASH:=Miniconda2-4.5.4-Linux-x86_64.sh
endif
CONDAURL:=https://repo.continuum.io/miniconda/$(CONDASH)
CONDADIR:=$(shell python -c 'import os; print(os.path.realpath("conda"))')
CONDA_ACTIVATE:=$(CONDADIR)/bin/activate
conda:
	wget "$(CONDAURL)" && \
	bash "$(CONDASH)" -b -p conda && \
	rm -f "$(CONDASH)" && \
	unset PYTHONHOME; unset PYTHONPATH && \
	source "$(CONDA_ACTIVATE)" && \
	conda config --add channels defaults && \
	conda config --add channels bioconda 

conda-install:
	unset PYTHONHOME; unset PYTHONPATH && \
	source "$(CONDA_ACTIVATE)" && \
	conda install -y r-base=3.4.3 \
    r=3.4.3 \
    r-ggplot2=2.2.1 \
    gdc-client=1.3.0
# r-data.table=1.10.4

conda-search:
	unset PYTHONHOME; unset PYTHONPATH && \
	source "$(CONDA_ACTIVATE)" && \
	conda search 'gdc-client'
# NOTE: methods of getting local conda into PATH;
# unset PYTHONHOME; unset PYTHONPATH; export PATH=$(CONDADIR)/bin:$$PATH; \
# unset PYTHONHOME; unset PYTHONPATH; source "$(CONDA_ACTIVATE)" && \
# TIMESTAMP=$(shell date +%s)
# clean-conda:
# [ -d conda ] && mv conda ".conda.old.${TIMESTAMP}" && rm -rf ".conda.old.${TIMESTAMP}" &

# ~~~~~ SETUP PIPELINE ~~~~~ #
./nextflow:
	export NXF_VER="$(NXF_VER)" && \
	curl -fsSL get.nextflow.io | bash

install: ./nextflow

update: ./nextflow
	./nextflow self-update

# ~~~~~ RUN PIPELINE ~~~~~ #
gdc: install
	./nextflow run gdc.nf -profile phoenix $(EP)

run: install
	./nextflow run main.nf -resume -profile phoenix $(EP)

# ~~~~~ CLEANUP ~~~~~ #
clean-traces:
	rm -f trace*.txt.*

clean-logs:
	rm -f .nextflow.log.*

clean-reports:
	rm -f *.html.*

clean-flowcharts:
	rm -f *.dot.*

clean-output:
	[ -d output ] && mv output oldoutput && rm -rf oldoutput &

clean-work:
	[ -d work ] && mv work oldwork && rm -rf oldwork &

# deletes files from previous runs of the pipeline, keeps current results
clean: clean-logs clean-traces clean-reports clean-flowcharts

# deletes all pipeline output
clean-all: clean clean-output clean-work
	[ -d .nextflow ] && mv .nextflow .nextflowold && rm -rf .nextflowold &
	rm -f .nextflow.log
	rm -f *.png
	rm -f trace*.txt*
	rm -f *.html*
