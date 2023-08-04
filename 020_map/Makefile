HISAT2 :=  docker run --rm -u $(id -u):$(id -g) -v $(PWD):/trans -w /trans/020_map \
  dceoy/hisat2

%.paired.sam: ../010_trim/%_1.paired.fastq.gz ../010_trim/%_2.paired.fastq.gz
  $(HISAT2) -q -p8 --dta -x ../000_input/genn.index/grch38_rep/genome_rep \
  -1 $(word 1,$^) \
  -2 $(word 2,$^) \
  -S $@
