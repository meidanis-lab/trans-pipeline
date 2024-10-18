TRIMMO :=  docker run --rm -u $(id -u):$(id -g) -v $(PWD):/trans -w /trans/010_trim chrishah/trimmomatic-docker trimmomatic

%_1.paired.fastq.gz: ../000_input/fastq/%_1.fastq.gz ../000_input/fastq/%_2.fastq.gz
	$(TRIMMO) PE $^ \
	$*_1.paired.fastq.gz $*_1.unpaired.fastq.gz \
	$*_2.paired.fastq.gz $*_2.unpaired.fastq.gz \
	ILLUMINACLIP:/usr/src/Trimmomatic/0.38/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10:2:TRUE \
	LEADING:33 TRAILING:33 MINLEN:36
