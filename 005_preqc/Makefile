FASTQC :=  docker run --rm -u $(id -u):$(id -g) -v $(PWD):/trans -w /trans/005_preqc \
	staphb/fastqc:0.11.9 fastqc

%_fastqc.html: ../000_input/fastq/%.fastq.gz
	$(FASTQC) $<  -o .
