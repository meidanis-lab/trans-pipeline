ASSEMBLE :=  docker run --rm -u $(id -u):$(id -g) -v $(PWD):/trans -w /trans/030_assemble \
	bschiffthaler/stringtie 

%.paired.gtf: ../025_sam/%.paired.bam ../000_input/ref.genome/Homo_sapiens.GRCh38.84.gtf
	$(ASSEMBLE) -o $@ -G $(word 2,$^) $<
