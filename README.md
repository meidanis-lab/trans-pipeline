# trans-pipeline
Pipeline for gene expression profile analisis of transcripts

Esse repositório visa compartilhar os códigos desenvolidos para análise de expressão diferencial de transcritos a partir de matrizez de contadores. 
Os dados de entrada, como os counts, metados e valores de fpkm, para os scripts podem ser baixados pelo repositório Zenodo (https://doi.org/10.5281/zenodo.7586587). 

O arquivo Dockerfile é um arquivo do Docker para criar uma imagem no RStudio com todos os pacotes necessários para rodar os scripts. 

O arquivo pipeline_deseq2 realiza análise de expressão diferencial entre grupos, constrói Diagrama de Venn e heatmap, para o programa DESeq2. 
O mesmo vale para o arquivo pipeline_edger, porém, para o programa edgeR. 

O arquivo pipeline_fpkm faz análises de variância entre os valores de FPKM de genes e suas amostras. 
