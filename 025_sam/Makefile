SAMTOOLS :=  docker run --rm -u $(id -u):$(id -g) -v $(PWD):/trans -w /trans/025_sam \
  staphb/samtools samtools

%.paired.bam: ../020_map/%.paired.sam
  $(SAMTOOLS) sort -o $@ $<
