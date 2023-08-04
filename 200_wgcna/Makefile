RSTUDIO_NET := docker run -v $(PWD)/:/trans -w /trans/200_wgcna --rm -e PASSWORD=212010 -p 80:8787 -e USERID=$UID trans:latest

nodes.csv: ../100_counts/counts.csv ../105_meta/meta.csv design.txt ../115_deg/deg_list
	$(RSTUDIO_NET) Rscript ./auto_net.R $(word 1,$^) $(word 2,$^) $(shell cat $(word 3,$^)) $(word 4,$^) $@
