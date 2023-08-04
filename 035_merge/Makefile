MERGE :=  docker run --rm -u $(id -u):$(id -g) -v $(PWD):/trans -w /trans/035_merge \
	bschiffthaler/stringtie

merge_list.gtf: merge_list.txt ../000_input/ref.genome/Homo_sapiens.GRCh38.84.gtf
	$(MERGE) --merge -G $(word 2,$^) -o $@ $<

merge_list.txt: $(wildcard ../030_assemble/*.paired.gtf)
	echo $^ | tr " " "\n" > $@
