SHELL:=/bin/bash
.ONESHELL:
ENV_SH:=env.sh
export PATH:=$(CURDIR):$(PATH)
export NXF_VER:=20.10.0
./nextflow:
	. "$(ENV_SH)" ; \
	curl -fsSL get.nextflow.io | bash

install: ./nextflow

# original container;
# docker://cmopipeline/getbasecountsmultisample:1.2.2

# use this one instead beacuse we need bash for nextflow
# https://hub.docker.com/layers/mskcc/helix_filters_01/getbasecountsmultisample-1.2.2/images/sha256-32c64ff26ce85c791cadb597a2b665ef516e396ace710ffe3cdc155b95c9d014?context=repo
# mskcc/helix_filters_01:getbasecountsmultisample-1.2.2
getbasecountsmultisample.1.2.2.sif:
	. "$(ENV_SH)" ; \
	singularity pull --force --name getbasecountsmultisample.1.2.2.sif docker://mskcc/helix_filters_01:getbasecountsmultisample-1.2.2

helix_filters_01.latest.sif:
	. "$(ENV_SH)" ; \
	singularity pull --force --name helix_filters_01.latest.sif docker://mskcc/helix_filters_01:latest

singularity-pull: getbasecountsmultisample.1.2.2.sif helix_filters_01.latest.sif

run:
	. "$(ENV_SH)" ; \
	nextflow run main.nf

bash:
	. "$(ENV_SH)" ; \
	bash

clean:
	rm -f .nextflow.log.*
	rm -f *.html.*
	rm -f trace.txt.*

clean-all: clean
	rm -f .nextflow.log*
	rm -f *.html
	rm -f trace.txt*
	rm -rf work
	rm -rf output
