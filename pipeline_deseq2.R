#Expressão diferencial dos subtipos
library("VennDiagram")
library("BiocManager")
library("ggplot2")
library("data.table")
library("dplyr")
library("tidyr")
library("plotly")
library("stringr")
library("DESeq2")
library("readr")
library("ComplexHeatmap")
library("circlize")

#Trazendo counts e metadata
### Passar nomes dos arquivos como parâmetro
counts <- read.csv("2022-transcriptome/diff-exp/input/counts.csv", sep = " ")
kidney.meta <- read.csv("2022-transcriptome/diff-exp/input/kidney.meta.csv", sep = " ")

#Resolvendo o problema ddo DESeq2: 
#Tirando os caracteres que o DESeq2 não aceita
kidney.meta <- as.data.frame(apply(kidney.meta, 
                                   2, function(kidney.meta) gsub("\\s", ".", 
                                                                 kidney.meta)))
kidney.meta$tumor_type <- as.character(kidney.meta$tumor_type)
kidney.meta$tumor_type[kidney.meta$tumor_type == "Metabolically.Divergent.(MD-)ChRCC"] <- "Metabolically.Divergent.ChRCC"
kidney.meta$tumor_type <- as.factor(kidney.meta$tumor_type)

#################################################
#Construindo o dataset
#Resolvendo o problema de linearidade das amostras 
kidney.meta$group <- paste0(gsub(" ", ".", 
                                 kidney.meta$tumor_type), 
                            "_", 
                            kidney.meta$pathology_kidney)

#Convertendo o design em fator
kidney.meta$group <- as.factor(kidney.meta$group)
kidney.meta$pathology_kidney <- as.factor(kidney.meta$pathology_kidney)
kidney.meta$tumor_type <- as.factor(kidney.meta$tumor_type)
dds.top <- DESeqDataSetFromMatrix(countData = counts, 
                                  colData = kidney.meta, 
                                  design = ~group)

#Mantendo aqueles com mais de 10 de contagem
keep <- rowSums(counts(dds.top)) >= 10
dds.top <- dds.top[keep,]

#Analise DE
DE.top <- DESeq(dds.top)
RE.top <- results(DE.top)

###CCRCC X T1
CC.T1 <- results(DE.top,
                 contrast = c("group", 
                              "Clear.cell.RCC_ccRCC",
                              "Type.1.Papillary.RCC_PRCC"))

CC.T1 <- subset(CC.T1, 
                CC.T1$padj < 0.05 & 
                  (CC.T1$baseMean > 50) & 
                  (abs(CC.T1$log2FoldChange) > 2))

CC.T1$sig <- as.factor(ifelse(CC.T1$padj < 0.05 & 
                                abs(CC.T1$log2FoldChange) >2 
                              & CC.T1$baseMean > 50, 
                              ifelse(CC.T1$log2FoldChange > 2, 
                                     'up','down'), 'not'))
summary(CC.T1$sig)
#write.table(CC.T1, file = "CC.T1.csv")

CC.T1.up <- rownames(CC.T1)[which(CC.T1$sig =="up")]
CC.T1.down <- rownames(CC.T1)[which(CC.T1$sig =="down")]


####CCRCC X CHRCC
CC.CHRCC <- results(DE.top,
                    contrast = c("group", 
                                 "Clear.cell.RCC_ccRCC",
                                 "Chromophobe.renal.cell.carcinoma_ChRCC"))

CC.CHRCC <- subset(CC.CHRCC, 
                   CC.CHRCC$padj < 0.05 & 
                     (CC.CHRCC$baseMean > 50) & 
                     (abs(CC.CHRCC$log2FoldChange) > 2))

CC.CHRCC$sig <- as.factor(ifelse(CC.CHRCC$padj < 0.05 & 
                                   abs(CC.CHRCC$log2FoldChange) >2 
                                 & CC.CHRCC$baseMean > 50, 
                                 ifelse(CC.CHRCC$log2FoldChange > 2, 
                                        'up','down'), 'not'))
summary(CC.CHRCC$sig)
#write.table(CC.CHRCC, file = "CC.CHRCC.csv")

CC.CHRCC.up <- rownames(CC.CHRCC)[which(CC.CHRCC$sig =="up")]
CC.CHRCC.down <- rownames(CC.CHRCC)[which(CC.CHRCC$sig =="down")]

####T1 X CHRCC
T1.CHRCC <- results(DE.top,
                    contrast = c("group", 
                                 "Type.1.Papillary.RCC_PRCC",
                                 "Chromophobe.renal.cell.carcinoma_ChRCC"))

T1.CHRCC <- subset(T1.CHRCC, 
                   T1.CHRCC$padj < 0.05 & 
                     (T1.CHRCC$baseMean > 50) & 
                     (abs(T1.CHRCC$log2FoldChange) > 2))

T1.CHRCC$sig <- as.factor(ifelse(T1.CHRCC$padj < 0.05 & 
                                   abs(T1.CHRCC$log2FoldChange) >2 
                                 & T1.CHRCC$baseMean > 50, 
                                 ifelse(T1.CHRCC$log2FoldChange > 2, 
                                        'up','down'), 'not'))
summary(T1.CHRCC$sig)
#write.table(T1.CHRCC, file = "T1.CHRCC.csv")

T1.CHRCC.up <- rownames(T1.CHRCC)[which(T1.CHRCC$sig =="up")]
T1.CHRCC.down <- rownames(T1.CHRCC)[which(T1.CHRCC$sig =="down")]

#####Diagrama de Venn
venn.up.top <- venn.diagram(x = list(CC.T1.up,
                                     CC.CHRCC.up,
                                     T1.CHRCC.up),
                            category.names = c("ccRCC vs PRCC 1",
                                               "ccRCC vs ChRCC",
                                               "PRCC 1 vs ChRCC"),
                            main = "Genes Up Regulated - DESeq2",main.just = c(2,2),
                            main.cex = 2,
                            fill = c("red", "green", "blue"),
                            alpha = c(0.5, 0.5,0.5), cex = 2,cat.fontface = 2,
                            lty =1, filename=NULL, cat.cex=2, 
                            cat.just=list(c(1,-2) , c(0.5,-34) , c(0.2,-2.6)))

grid.newpage()
grid.draw(venn.up.top)

venn.down.top <- venn.diagram(x = list(CC.T1.down,
                                       CC.CHRCC.down,
                                       T1.CHRCC.down), 
                              category.names = c("ccRCC vs PRCC 1",
                                                 "ccRCC vs ChRCC",
                                                 "PRCC 1 vs ChRCC"),
                              main = "Genes Down Regulated - DESeq2", main.just = c(2,2),
                              main.cex = 2,
                              fill = c("red", "green", "blue"),
                              alpha = c(0.5, 0.5,0.5), cex = 2,cat.fontface = 2,
                              lty =1, filename=NULL, cat.cex=2, 
                              cat.just=list(c(1,-2) , c(0.5,-34) , c(0.2,-2.6)))
grid.newpage()
grid.draw(venn.down.top)

##############################


###DOWN genes
top5.down <- Reduce(intersect, list(CC.T1.down, CC.CHRCC.down, T1.CHRCC.down))

top520.down <- Reduce(intersect, list(T1.CHRCC.down, CC.CHRCC.down))

top22.down <- Reduce(intersect, list(CC.CHRCC.down, CC.T1.down))

top537.down <- union(top22.down, top520.down)

#### Up genes
top8.up <- Reduce(intersect, list(CC.T1.up, CC.CHRCC.up, T1.CHRCC.up))

top331.up <- Reduce(intersect, list(T1.CHRCC.up, CC.CHRCC.up))

top73.up <- Reduce(intersect, list(CC.CHRCC.up, CC.T1.up))

top396.up <- union(top73.up, top331.up)

#########Juntando
deg <- c(top537.down, top396.up)

###Normalização vst

vsd <- vst(DE.top, blind = F)
assay.deg <- assay(vsd)[deg,]

coldata <- colData(vsd)[,c("pathology_kidney", "tumor_type")]
rows.T1 <- rownames(coldata)[which(coldata$tumor_type=="Type.1.Papillary.RCC")]
rows.CC <- rownames(coldata)[which(coldata$tumor_type == "Clear.cell.RCC")]
rows.CH <- rownames(coldata)[which(coldata$tumor_type == "Chromophobe.renal.cell.carcinoma")]
samples.deg <- c(rows.T1, rows.CC, rows.CH)
coldata <- coldata[samples.deg,]

assay.deg <- assay.deg[,samples.deg]
min(assay.deg)
max(assay.deg)

all(colnames(assay.deg) == rownames(coldata))

####
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun(seq(-3, 3))

###Anotação
type <- coldata$tumor_type
top.anno <- HeatmapAnnotation(foo = as.matrix(type),
                              col = list(foo = c("Clear.cell.RCC" =  "light green", 
                                                "Chromophobe.renal.cell.carcinoma" = "violet",
                                                "Type.1.Papillary.RCC" = "dark gray")),
                              annotation_label = c("Cancer subtype"),
                              annotation_name_gp = gpar(fontsize = 15),
                              annotation_legend_param = list(labels = c("ChRCC - Unclassified",
                                                                        "ccRCC",
                                                                        "PRCC - Type 1")))
###HeatMap
Heatmap(assay.deg, 
        show_row_names = F,
        show_column_names = F,
        column_split = coldata$tumor_type, 
        top_annotation = top.anno,
        row_title = "Genes", row_title_gp = gpar(fontsize = 10),
        column_title = "Heatmap VST (variance stabilizig transformation)", 
        column_title_gp = gpar(fontsize = 16), 
        heatmap_legend_param = list(col_fun = col_fun, 
                                    title = "VST", 
                                    at = c(4,10,20)))




