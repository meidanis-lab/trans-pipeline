#!/bin/bash

#100_counts
make -C 100_counts/ counts.csv

#105_meta
make -C 105_meta/ contrast_test.txt
rm 105_meta/output_file.tmp 105_meta/temp_file.txt 105_meta/output.txt

#110_contrast
make -C 110_contrast/ output

#115_deg
make -C 115_deg/ deg_list
