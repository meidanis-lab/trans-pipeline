SHELL := /bin/bash

RSTUDIO_DATASEQ := docker run -v $(PWD)/..:/trans -w /trans/110_contrast --rm -e PASSWORD=212010 -p 80:8787 -e USERID=$UID trans:latest

output: ../105_meta/contrasts.txt
	for colctrcnd in $(shell cat $<); do \
		touch $$colctrcnd; make output_$$colctrcnd; \
		rm $$colctrcnd; \
	done

output_%: ../100_counts/counts.csv ../105_meta/meta.csv %
	echo "Processing $*"
	$(RSTUDIO_DATASEQ) Rscript ./de_analysis.R $(word 1,$^) $(word 2,$^) `echo $* | cut -d ',' -f1` `echo $* | cut -d ','  -f2` `echo $* | cut -d ',' -f3` $@

