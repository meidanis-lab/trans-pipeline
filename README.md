# trans-pipeline
Pipeline for gene expression profile analisis of transcripts

This repository aims to share the codes developed for the analysis of differential expression of transcripts using data from RNA-seq. The input data, such as counts, metadata, and fpkm values, can be downloaded from the Zenodo repository (https://doi.org/10.5281/zenodo.7586587).

The Dockerfile present in the docs directory is a Docker file used to create an RStudio image with all the necessary packages to run the scripts.

The `parcial_report` directory includes the scripts written for the partial project report. It complements the `pipeline_deseq2` file, which performs differential expression analysis between groups, constructs Venn diagrams and heatmaps for the DESeq2 program. The same applies to the `pipeline_edger` file, but for the edgeR program. The `pipeline_fpkm` file performs variance analyses between FPKM values of genes and their samples related to renal carcinoma data.

The remaining directories (`005_preqc` to `205_gwena`) contain the necessary scripts for sample pre-processing, differential expression analysis, network analysis, and functional annotation. For a better understanding, you can access the tutorial at https://rpubs.com/jpietro/vignette_pipeline.

