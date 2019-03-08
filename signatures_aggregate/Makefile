SHELL:=/bin/bash

# no default action to take
none:

# run this first; make a local copy of all the files we need from the pipelines
copy:
	find /gpfs/data/molecpathlab/production/NGS580/181228_NB501073_0088_AH5NJWBGX9 -type f -name "*signatures.Rds" -path "*/output/*" -exec cp -va {} data/181228_NB501073_0088_AH5NJWBGX9/output/signatures/ \;
	find /gpfs/data/molecpathlab/production/NGS580/190201_NB501073_0092_AH7333BGX9 -type f -name "*signatures.Rds" -path "*/output/*" -exec cp -va {} data/190201_NB501073_0092_AH7333BGX9/output/signatures/ \;
	find /gpfs/data/molecpathlab/production/NGS580/190208_NB501073_0094_AHFHK3BGX9 -type f -name "*signatures.Rds" -path "*/output/*" -exec cp -va {} data/190208_NB501073_0094_AHFHK3BGX9/output/signatures/ \;

# scripts to run
WEIGHTS_SCRIPT:=/gpfs/data/molecpathlab/projects/grab_signatures_quick_for_George/weights.R
PASTE_SCRIPT:=/gpfs/data/molecpathlab/projects/grab_signatures_quick_for_George/paste-col.py
CONCAT_SCRIPT:=/gpfs/data/molecpathlab/projects/grab_signatures_quick_for_George/concat-tables.py
PWD:="$(shell pwd)"

# Singularity container that contains R and deconstructSigs
CONTAINER:=/gpfs/data/molecpathlab/containers/NGS580-nf/deconstructSigs-1.8.0.simg

# get the list of input files to process
INPUT_FILES:=$(shell find data/ -type f -name "*_signatures.Rds")
.PHONY: $(INPUT_FILES)
run: $(INPUT_FILES)

# processes each file in the list
# load Singularity module
# enter Singularity container with directories mounted, run commands to extract weights.tsv
# add a Sample ID column to the tsv
$(INPUT_FILES):
	input_file="$$(readlink -f "$@")" ; \
	input_base="$$(basename $${input_file})" ; \
	tmp_file="$$(echo "$@" | sed -e 's|^data/|output/|g' -e 's|.Rds$$|.weights.tmp|g' )" ; \
	tmp_file="$$(readlink -m $${tmp_file} )" ; \
	output_file="$$(echo "$@" | sed -e 's|^data/|output/|g' -e 's|.Rds$$|.weights.tsv|g' )" ; \
	output_file="$$(readlink -m $${output_file} )" ; \
	mkdir -p "$$(dirname $${output_file:?})" ; \
	module load singularity/2.5.2 ; \
	echo "[files] $${input_file:?} , $${tmp_file:?} , $${output_file:?}" ; \
	singularity exec -B $${PWD} \
	-B $$(dirname $${input_file}) \
	-B $$(readlink -f data) \
	-B $$(readlink -f output) \
	$(CONTAINER) \
	/bin/bash -c "cd $(PWD); $(WEIGHTS_SCRIPT) $${input_file} $${tmp_file} " && \
	$(PASTE_SCRIPT) -i "$${tmp_file}" -o "$${output_file}" --header "File" -v "$${input_base}"

# find the TSV files and concat them all with the concat script
TSV_FILES:=$(shell find output/ -type f -name "*_signatures.weights.tsv")
concat:
	$(CONCAT_SCRIPT) $(TSV_FILES) > weights.tsv
