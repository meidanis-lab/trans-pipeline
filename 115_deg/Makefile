deg_list: $(wildcard ../140_contrast/output_*)*
	find $^  -type f -name "output_*" -exec sed  -E '/^(not|down|up)\s/ s/^(\S+).*/\1/p' {} \; \
		| awk '{print $$1,$$8}' | sort | uniq > $@
