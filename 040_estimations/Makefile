ESTIMATION  :=  docker run --rm -u $(id -u):$(id -g) -v $(PWD):/trans -w /trans/040_estimations \
	bschiffthaler/stringtie

%.paired.gtf: ../025_sam/%.paired.bam ../035_merge/merge_list.gtf
	$(ESTIMATION) -e -B -G $(word 2,$^) -o $@ $<
