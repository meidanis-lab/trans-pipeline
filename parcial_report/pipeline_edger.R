library("ggplot2")
library("data.table")
library("stringr")
library("readr")
library("pheatmap")
library("ComplexHeatmap")
library("magick")
library("edgeR")
library("VennDiagram")
library("BiocManager")
library("dplyr")
library("tidyr")
library("plotly")
library("circlize")

#######PIPELINE ARTIGO#################################
#Reading counts
counts <- read.csv("2022-transcriptome/diff-exp/input/counts.csv", sep = " ")

#Metadados só para o edgeR
metadata <- read.csv("2022-transcriptome/diff-exp/input/kidney.meta.csv", sep = " ")
metadata <- as.data.frame(apply(metadata, 
                                2, function(metadata) gsub("\\s", ".", 
                                                           metadata)))
metadata$tumor_type <- as.character(metadata$tumor_type)
metadata$tumor_type[metadata$tumor_type == "Metabolically.Divergent.(MD-)ChRCC"] <- "Metabolically.Divergent.ChRCC"
metadata$tumor_type <- as.factor(metadata$tumor_type)

#design
group <- factor(paste(metadata$pathology_kidney, 
                      metadata$tumor_type, sep = "_"))
design <- model.matrix(~0+group)
colnames(design) <- levels(group)

#DGE list
DGElist <- DGEList(counts, group = group)

#DGE filter
keep <- rowSums(cpm(DGElist) > 100) >= 2
table(keep)
DGElist <- DGElist[keep, , keep.lib.sizes = F]

#Normalização
DGEnorm <- calcNormFactors(DGElist)
DGEnorm$samples

#Dispersão
y <- estimateDisp(DGEnorm, design, robust = T)
fit <- glmQLFit(y, design, robust=T)
summary(fit$df.prior)


## GRUPOS COM MAIS AMOSTRAS SEM CONSIDERAR CLASSIFICAÇÃO
###############CCRCC X PRCC T1#######################
CC.T1 <- makeContrasts(ccRCC_Clear.cell.RCC - PRCC_Type.1.Papillary.RCC,
                         levels = design)

FT.CC.T1 <- glmQLFTest(fit, contrast = CC.T1)
topTags(FT.CC.T1)
summary(decideTests(FT.CC.T1))

TR.CC.T1 <- glmTreat(fit, contrast=CC.T1, lfc=log2(1.5))
summary(decideTests(TR.CC.T1))

CC.T1 <- TR.CC.T1$table
CC.T1<- CC.T1[abs(CC.T1$logFC) > 2,]
CC.T1$sigs <- as.factor(
  ifelse(abs(CC.T1$logFC) > 2,
         ifelse(CC.T1$logFC > 2 ,'up','down'),'not'))
summary(CC.T1$sigs)
#write.table(CC.T1, file = "CC.T1.csv")

CC.T1.up <- rownames(CC.T1)[which(CC.T1$sig =="up")]
CC.T1.down <- rownames(CC.T1)[which(CC.T1$sig =="down")]

###############CCRCC X CHRCC#######################
CC.CH <- makeContrasts(ccRCC_Clear.cell.RCC - ChRCC_Chromophobe.renal.cell.carcinoma,
                       levels = design)

FT.CC.CH <- glmQLFTest(fit, contrast = CC.CH)
topTags(FT.CC.CH)
summary(decideTests(FT.CC.CH))

TR.CC.CH <- glmTreat(fit, contrast=CC.CH, lfc=log2(1.5))
summary(decideTests(TR.CC.CH))

CC.CH <- TR.CC.CH$table
CC.CH<- CC.CH[abs(CC.CH$logFC) > 2,]
CC.CH$sigs <- as.factor(
  ifelse(abs(CC.CH$logFC) > 2,
         ifelse(CC.CH$logFC > 2 ,'up','down'),'not'))
summary(CC.CH$sigs)
#write.table(CC.T1, file = "CC.CH.csv")

CC.CH.up <- rownames(CC.CH)[which(CC.CH$sig =="up")]
CC.CH.down <- rownames(CC.CH)[which(CC.CH$sig =="down")]

###############CCRCC X CHRCC#######################
T1.CH <- makeContrasts(PRCC_Type.1.Papillary.RCC - ChRCC_Chromophobe.renal.cell.carcinoma,
                       levels = design)

FT.T1.CH <- glmQLFTest(fit, contrast = T1.CH)
topTags(FT.T1.CH)
summary(decideTests(FT.T1.CH))

TR.T1.CH <- glmTreat(fit, contrast=T1.CH, lfc=log2(1.5))
summary(decideTests(TR.T1.CH))

T1.CH <- TR.T1.CH$table
T1.CH <- T1.CH[abs(T1.CH$logFC) > 2,]
T1.CH$sigs <- as.factor(
  ifelse(abs(T1.CH$logFC) > 2,
         ifelse(T1.CH$logFC > 2 ,'up','down'),'not'))
summary(T1.CH$sigs)
#write.table(CC.T1, file = "T1.CH")

T1.CH.up <- rownames(T1.CH)[which(T1.CH$sig =="up")]
T1.CH.down <- rownames(T1.CH)[which(T1.CH$sig =="down")]

#####Diagrama de Venn
venn.up.top <- venn.diagram(x = list(CC.T1.up,
                                     CC.CH.up,
                                     T1.CH.up), 
                            category.names = c("ccRCC vs PRCC 1",
                                               "ccRCC vs ChRCC",
                                               "PRCC 1 vs ChRCC"),
                            main = "Genes Up Regulated - edgeR",main.just = c(2,2),
                            main.cex = 2,
                            fill = c("red", "green", "blue"),
                            alpha = c(0.5, 0.5,0.5), cex = 2,cat.fontface = 2,
                            lty =1, filename=NULL, cat.cex=2, 
                            cat.just=list(c(1,-2) , c(0.5,-34) , c(0.2,-2.6)))
grid.newpage()
grid.draw(venn.up.top)

venn.down.top <- venn.diagram(x = list(CC.T1.down,
                                       CC.CH.down,
                                       T1.CH.down), 
                              category.names = c("ccRCC vs PRCC 1",
                                                 "ccRCC vs ChRCC",
                                                 "PRCC 1 vs ChRCC"),
                              main = "Genes Down Regulated - edgeR", main.just = c(2,2),
                              main.cex = 2,
                              fill = c("red", "green", "blue"),
                              alpha = c(0.5, 0.5,0.5), cex = 2,cat.fontface = 2,
                              lty =1, filename=NULL, cat.cex=2, 
                              cat.just=list(c(1,-2) , c(0.5,-34) , c(0.2,-2.6)))
grid.newpage()
grid.draw(venn.down.top)


#####GENES UP
top10.up <- Reduce(intersect, list(CC.T1.up, CC.CH.up, T1.CH.up))

top76.up <- Reduce(intersect, list(CC.CH.up, CC.T1.up))
top231.up <- Reduce(intersect, list(T1.CH.up,CC.CH.up))

top297.up <- union(top76.up, top231.up)

###### GENES DOWN
top5.down <- Reduce(intersect, list(CC.T1.down, CC.CH.down, T1.CH.down))
top18.down <- Reduce(intersect, list(CC.CH.down, CC.T1.down))
top402.down <- Reduce(intersect, list(T1.CH.down, CC.CH.down))

top415.down <- union(top18.down,top402.down)

##### Juntando
deg.edger <- c(top297.up, top415.down)

###Selecionando as amostras dos DEG
rows.T1 <- rownames(metadata)[which(metadata$tumor_type=="Type.1.Papillary.RCC")]
rows.CC <- rownames(metadata)[which(metadata$tumor_type == "Clear.cell.RCC")]
rows.CH <- rownames(metadata)[which(metadata$tumor_type == "Chromophobe.renal.cell.carcinoma")]
samples.deg <- c(rows.T1, rows.CC, rows.CH)
metadata <- metadata[samples.deg,]

logcpm <- cpm(y, log = T)
logcpm <- logcpm[deg.edger,samples.deg]

all(rownames(metadata) == colnames(logcpm))

#Anotação
type.edger <- metadata$tumor_type
anno <- HeatmapAnnotation(foo = as.matrix(type.edger),
                              col = list(foo = c("Clear.cell.RCC" =  "light green", 
                                                 "Chromophobe.renal.cell.carcinoma" = "violet",
                                                 "Type.1.Papillary.RCC" = "dark gray")),
                              annotation_label = c("Cancer subtype"),
                              annotation_name_gp = gpar(fontsize = 15),
                          annotation_legend_param = list(labels = c("ChRCC - Unclassified",
                                                                    "ccRCC",
                                                                    "PRCC - Type 1")))

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun(seq(-3, 3))

###HeatMap
min(logcpm)
max(logcpm)

Heatmap(logcpm, 
        show_row_names = F,
        show_column_names = F,
        column_split = metadata$tumor_type,
        row_title = "Genes", row_title_gp = gpar(fontsize = 15),
        column_title = "Heatmap logCPM", column_title_gp = gpar(fontsize = 16),
        top_annotation = anno,
        heatmap_legend_param = list(col_fun = col_fun, 
                                    title = "log CPM"))


