POSFASTQC :=  docker run --rm -u $(id -u):$(id -g) -v $(PWD):/trans -w /trans/015_posqc \
  staphb/fastqc:0.11.9 fastqc

%.paired.fastq.html: ../010_trim/%.paired.fastq.gz
  $(POSFASTQC) $< -o .
