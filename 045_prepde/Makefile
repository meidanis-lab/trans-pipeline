count_matrix.csv: ID
	./prepDE.py3 -i $< -t $@

ID: names.txt paths.txt
	paste $(word 1,$^) $(word 2,$^) > $@

paths.txt: $(wildcard ../040_estimations/*.paired.gtf)
	echo $^ | tr " " "\n" > $@

names.txt: paths.txt
	sed 's|.*/||;s/\.paired.gtf//' $< > $@
