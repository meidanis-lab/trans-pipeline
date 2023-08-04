args = commandArgs(trailingOnly=TRUE)

### Usage:
### Rscript --vanilla de_analysis.R <counts.csv> <meta.csv> <column> <control> <condition> <output>

countsfile = args[1]
metafile = args[2]
column = args[3]
print(column)
degfile = args[4]
nodesfiles = args[5]

library(WGCNA)
library(DESeq2)
library(enrichplot)
library(biomaRt)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(clusterProfiler)
library(AnnotationDbi)


blockwiseModules
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

###############DATA INPUT, CLEANING AND NORMALIZATION####################

#Reading the files
counts = read.csv(countsfile, header = T, row.names = 1, sep = ",")
meta = read.csv(metafile, header = T, row.names = 1, sep = ",")

#Checking if row names matches column names 
all(rownames(meta)==colnames(counts))

#Building the dataset from matrix
dds = DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = as.formula(paste("~",column)))
keep <-  rowSums(counts(dds)) >= 100
dds <- dds[keep,]
vsd <- vst(dds, blind = T)
datExpr0 = assay(vsd)

png(file = "plotpca.png", width = 800, height = 400)
plotPCA(vsd, intgroup= colnames(meta))
dev.off()

#Turning categorical metadata to numerical 
samples = rownames(meta)
meta <- sapply(meta, function(x) as.numeric(factor(x)))
rownames(meta) = samples

# Calculate the variance or standard deviation for each gene
gene_variability <- apply(datExpr0, 1, var)  
# Determine the number of genes to select (5% of total genes)
num_genes <- round(0.05 * length(gene_variability))
# Sort the genes based on their variability in descending order
sorted_genes <- order(gene_variability, decreasing = TRUE)
# Select the top 5% of genes with the highest variability
selected_genes <- sorted_genes[1:num_genes]
#Returning to datExpr0
datExpr0 = datExpr0[selected_genes,]
datExpr0 = t(datExpr0)

#Checking for good genes
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


##Detecting outliers
sampleTree = hclust(dist(datExpr0), method = "average");
cluster_heights = sampleTree[["height"]]
summary_clus = summary(sampleTree$height)

#Saving plot to detect outliers
sizeGrWindow(12,9)
png(file = "outliers.png", width = 800, height = 400)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()


# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = max(cluster_heights), minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
meta = meta[rownames(datExpr),]
rownames(meta)==rownames(datExpr)
sampleTree2 = hclust(dist(datExpr), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(meta, signed = FALSE)

#Saving plot of sample dendogram and trait heatmap
png(file = "dendogram.png", width = 800, height = 400)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(meta),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

######################BUILDING NETWORK#######################

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sft$powerEstimate

# Plot the results:
sizeGrWindow(9, 5)
png(file = "soft_power.png", width = 800, height = 400)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#Plotting net
cor <- WGCNA::cor
net = blockwiseModules(datExpr, power = sft$powerEstimate,
                       TOMType = "unsigned", minModuleSize = 150,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TOM",
                       verbose = 3)

cor <- stats::cor

##Detecting modules
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

####################RELATING WITH CLINICAL TRAITS######################

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, meta, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
png("heatmap.png", width = 800, height = 400)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
heatmap =labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(meta),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
print(heatmap)
dev.off()

module_MAX = rownames(moduleTraitCor)[which(moduleTraitCor == max(moduleTraitCor), 
                                            arr.ind = TRUE)[ ,1]]
module_MAX
module_MAX = substring(module_MAX,3)

##############ANNOTATION PART. 1#########################

#####Annotation with Biomart
ensembl_id = as.data.frame(colnames(datExpr))

listEnsembl()
ensembl = useEnsembl(biomart = "genes")
datasets = listDatasets(ensembl)

ensembl_con = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

attr = listAttributes(ensembl_con)
filters = listFilters(ensembl_con)

annotation = getBM(attributes = c("ensembl_transcript_id",
                                  "entrezgene_id",
                                  "external_gene_name",
                                  "ensembl_gene_id",
                                  "transcript_gencode_basic"), 
                   filters = "ensembl_transcript_id",
                   values = ensembl_id$`colnames(datExpr)`,
                   mart = ensembl_con)
write.table(annotation, file = "annotation.txt", sep = ",")

#####################TOM###############################
TOM = TOMsimilarityFromExpr(datExpr, power = 9);
# Select modules
module_t = c(module_MAX)
# Select module probes
probes = colnames(datExpr)
inModule = is.finite(match(moduleColors, module_t));
modProbes = probes[inModule];
modGenes = annotation$external_gene_name[match(modProbes, annotation$ensembl_transcript_id)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
modTOM <- as.data.frame(modTOM)

###Reading deg_list
deg_list = read.csv(degfile, sep = " ", header = T)
colnames(deg_list) = c("id", "sig")
deg.up <- deg_list$id[which(deg_list$sig=="up")]
deg.down <- deg_list$id[which(deg_list$sig=="down")]
deg_list = c(deg.up, deg.down)
deg_list = as.data.frame(unique(deg_list))
colnames(deg_list) = "id"


valid_ids = deg_list$id[deg_list$id %in% colnames(modTOM)]

modTOM_filter = modTOM[valid_ids, valid_ids]

##########NET WITH DEG
cyt = exportNetworkToCytoscape(modTOM_filter,
                               edgeFile = paste("edges-", paste(module_t, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("nodes-", paste(module_t, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.1, #0,05 0,15
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])
edges <- cyt$edgeData
colnames(edges)[1] <- "Source"
colnames(edges)[2] <- "Target"
edges<-edges[, -4:-6]

nodes <- cyt$nodeData
colnames(nodes)[1] <- "id"
colnames(nodes)[2] <- "Label"

hubs_modules = chooseTopHubInEachModule(datExpr,moduleColors,power = 9,type = "signed")

n_occur <- data.frame(table(nodes$Label))
n_occur[n_occur$Freq > 1,]

write.csv(edges,file = "edges.csv", row.names = F)
write.csv(nodes,file = nodesfiles, row.names = F)

##############GO E KEGG################

############ GO ######################
# Select modules
modules = c(module_MAX);
valid_entrez = annotation$entrezgene_id[match(valid_ids, annotation$ensembl_transcript_id)]
valid_ids
GO_results <- enrichGO(valid_entrez, "org.Hs.eg.db", ont = "BP")
fit = plot(barplot(GO_results, showCategory = 10, main = "GO"))
png("out_GO.png", res = 250,width = 1400, height = 2000)
print(fit)
dev.off()

############# KEGGREST ###############
valid_symbols = annotation$entrezgene_id[match(valid_ids, annotation$ensembl_transcript_id)]
keg = enrichKEGG(valid_symbols, keyType = "kegg", pvalueCutoff=0.5)
print(keg)
kegg = plot(barplot(keg, showCategory = 10))
png("kegg.png", res = 250, width = 1400, height = 2000)
print(kegg)
dev.off()

