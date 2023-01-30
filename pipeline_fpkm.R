library("BiocManager")
library("ggplot2")
library("data.table")
library("dplyr")
library("tidyr")
library("shiny")
library("plotly")
library("knitr")
library("stringr")
library("readr")
library("ComplexHeatmap")
library("circlize")
library("edgeR")
library("DESeq2")

###########FPKM E VARIÂNCIA
fpkm <- read.csv("2022-transcriptome/diff-exp/input/fpkm.csv", sep = " ") 

fpkm <- as.matrix(fpkm)

fpkm.order <- head(order(rowVars(fpkm), decreasing = T),500)

fpkm.var <- fpkm[fpkm.order,]

fpkm.log <- (log(fpkm.var+1))/0.12


#Metadados
kidney.meta <- read.csv("2022-transcriptome/diff-exp/input/kidney.meta.csv", sep = " ")

kidney.meta <- as.data.frame(apply(kidney.meta, 
                                   2, function(kidney.meta) gsub("\\s", ".", 
                                                                 kidney.meta)))
kidney.meta$tumor_type <- as.character(kidney.meta$tumor_type)
kidney.meta$tumor_type[kidney.meta$tumor_type == "Metabolically.Divergent.(MD-)ChRCC"] <- "Metabolically.Divergent.ChRCC"
kidney.meta$tumor_type <- as.factor(kidney.meta$tumor_type)

type <- kidney.meta$tumor_type

pathology <- kidney.meta$pathology_kidney

all(rownames(pathology) == colnames(fpkm.var))

#Anotação
anno_pathology <- HeatmapAnnotation(bar = as.matrix(pathology),
                             col = list(bar = c("PRCC" = "green",
                                                "ccRCC" = "red",
                                                "ChRCC" = "blue")),
                             annotation_label = c("Patologias"),
                             annotation_name_gp = gpar(fontsize = 10))

anno_type <- HeatmapAnnotation(foo = as.matrix(type),
                               col = list(foo =c("Clear.cell.RCC" =  "light green", 
                                                 "Chromophobe.renal.cell.carcinoma" = "violet",
                                                 "KIRP.CIMP" = "light blue", 
                                                 "Metabolically.Divergent.ChRCC" = "yellow",
                                                 "Type.1.Papillary.RCC" = "dark gray",
                                                 "Type.2.Papillary.RCC" = "brown",
                                                 "Unclassified.Papillary.RCC" = "red")),
                               annotation_label = c("Cancer subtype"),
                               annotation_name_gp = gpar(fontsize = 15),
                               annotation_legend_param = list(labels = c("ChRCC - Unclassified",
                                                                         "ccRCC",
                                                                         "PRCC - KIRP CIMP",
                                                                         "ChRCC - Metabolically Divergent",
                                                                         "PRCC - Type 1",
                                                                         "PRCC - Type 2",
                                                                         "PRCC - Unclassified")))

#Deixando o gráfico bonito
col_fun = colorRamp2(c(0, 50, 90), c("blue", "white", "red"))
col_fun(seq(-3, 3))

#Plotando
par(mfrow=c(2,2))
Heatmap(fpkm.log,
        show_row_names = F, 
        show_column_names = F,
        cluster_columns = T,
        cluster_rows = T,
        row_title = "Genes", row_title_gp = gpar(fontsize = 15),
        column_title = "Heatmap FPKM upper quartile normalized", column_title_gp = gpar(fontsize = 16), 
        top_annotation = anno_type,
        column_split = kidney.meta$tumor_type,
        heatmap_legend_param = list(col_fun = col_fun, 
                                    colorbar = "discrete",
                                    title = "Percentile", 
                                    at = c(0,40,90)))

