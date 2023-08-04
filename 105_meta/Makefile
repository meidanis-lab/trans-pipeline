SHELL := /bin/bash

contrast_test.txt: output_file.tmp
	cat $< | tr -d '"' > $@

output_file.tmp: output.txt temp_file.txt
	paste -d ',' "$(word 1,$^)" "$(word 2,$^)" > $@


temp_file.txt: ../105_meta/meta.csv
	header=$$(head -n 1 $<); \
	IFS=',' read -ra columns <<< "$$header"; \
	for ((i=1; i<$${#columns[@]}; i++)); do \
		echo -n "Column $$i ($${columns[i]}):"; \
		tail -n +2 "$<" | cut -d ',' -f $$((i+1)) | sort | uniq | paste -sd ',' - | sed 's/,$$//' >> "$@"; \
		echo; \
	done

output.txt: ../105_meta/meta.csv
	awk -F',' 'NR==1 {for(i=1;i<=NF;i++) print $$i}' $< | tr -d '"' | tail -n +2 >> $@
