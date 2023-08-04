library(WGCNA)
library(DESeq2)
library(enrichplot)
library(biomaRt)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

blockwiseModules
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
###############DATA INPUT, CLEANING AND NORMALIZATION####################

#Reading counts
counts = read.csv("100_counts/counts.csv", header = T, row.names = 1, sep = ",")
meta = read.csv("105_meta/meta.csv", header=T, row.names=1)

###Normalizing
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~Status)

keep <-  rowSums(counts(dds)) >= 100
dds <- dds[keep,]
vsd <- vst(dds, blind = T)
datExpr0 = assay(vsd)
pca = plotPCA(vsd, intgroup=c("Metastasis", "Status", "TP53", "CDK", "RB1"))

#My metada isn't numerical, so i need to convert into numerical values
# Create a mapping from categorical values to numerical values
mapping <- list(Status = c("Alive" = 0, "Dead" = 1),
                TP53 = c("no" = 0, "Mut" = 1),
                RB1 = c("no" = 0, "Mut" = 1),
                CDK = c("no" = 0, "Mut" = 1),
                COPS3 = c("no" = 0, "COP" = 1),
                Metastasis = c("no"=0, "Lung"=1, "Lung_bones"=2),
                PI3K.mTOR = c("no" = 0, "PTEN_Mut" = 1, "PTEN_del" = 2,
                              "NF1_Mut" = 3, "NF1_del" = 4, "AKT1_Mut" = 5,
                              "NF1_Mut.PDPK1_Mut" = 6))  
# Define a function to convert categorical values to numerical values
for (col in colnames(meta)) {
  if (col %in% names(mapping)) {
    meta[[col]] <- mapping[[col]][meta[[col]]]
  }
}


# Calculate the variance or standard deviation for each gene
gene_variability <- apply(datExpr0, 1, var)  # Using variance

# Determine the number of genes to select (4% of total genes)
num_genes <- round(0.05 * length(gene_variability))

# Sort the genes based on their variability in descending order
sorted_genes <- order(gene_variability, decreasing = TRUE)

# Select the top 4% of genes with the highest variability
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
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 200, col = "red")

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 800,minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
rownames(meta)==rownames(datExpr)

sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(meta, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(meta),
                    main = "Sample dendrogram and trait heatmap")


#save(datExpr, meta, file = "160_wgcna/DataInput.RData")

########################BUILDING NETWORK#######################

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sft$powerEstimate
#save(sft, file="160_wgcna/sft.RData")

# Plot the results:
sizeGrWindow(9, 5)
png(file = "plot_output.png", width = 800, height = 400)
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
net = blockwiseModules(datExpr, power =9,
                       TOMType = "unsigned", minModuleSize = 150,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "160_wgcna/TOM",
                       verbose = 3)

cor <- stats::cor

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
#save(MEs, moduleLabels, moduleColors, geneTree,
#     file = "NetworkConstruction-auto.RData")

###########RELATING WITH CLINICAL TRAITS#######
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, meta, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
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


##############################################
##Getting the names from each module
ids_blue = colnames(datExpr)[moduleColors=="blue"]
writeLines(ids_blue, "160_wgcna/ids_blue.txt")

ids_turquoise = colnames(datExpr)[moduleColors=="turquoise"]
writeLines(ids_turquoise, "160_wgcna/ids_turquoise.txt")

ids_brown = colnames(datExpr)[moduleColors=="brown"]
writeLines(ids_brown, "160_wgcna/ids_brown.txt")

ids_grey = colnames(datExpr)[moduleColors=="grey"]
writeLines(ids_grey, "160_wgcna/ids_grey.txt")

###################
deg_tp53 = read.csv("140_contrast/output_TP53,no,Mut", header = T, sep = "")
deg_rb1 = read.csv("140_contrast/output_RB1,no,Mut", header = T, sep = "")
deg_cdk = read.csv("140_contrast/output_CDK,no,Mut", header = T, sep = "")
writeLines(rownames(deg_tp53), "160_wgcna/deg_tp53")
writeLines(rownames(deg_rb1), "160_wgcna/deg_rb1")
writeLines(rownames(deg_cdk), "160_wgcna/deg_cdk")

deg_tp53.up <- rownames(deg_tp53)[which(deg_tp53$sig =="up")]
deg_tp53.down <- rownames(deg_tp53)[which(deg_tp53$sig =="down")]

deg_rb1.up <- rownames(deg_rb1)[which(deg_rb1$sig =="up")]
deg_rb1.down <- rownames(deg_rb1)[which(deg_rb1$sig =="down")]

deg_cdk.up <- rownames(deg_cdk)[which(deg_cdk$sig =="up")]
deg_cdk.down <- rownames(deg_cdk)[which(deg_cdk$sig =="down")]

intersect(deg_tp53.up, intersect(deg_cdk.up, deg_rb1.up))

deg_list = c(deg_tp53.up, deg_tp53.down, deg_rb1.up, deg_rb1.down, deg_cdk.up, deg_cdk.down)
deg_list = as.data.frame(unique(deg_list))
colnames(deg_list) = "id"


##########################################

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
write.table(annotation, file = "160_wgcna/annotation.txt", sep = ",")

# Recalculate topological overlap if needed
intersect(deg_list$id,ids_blue) #143
intersect(deg_list$id,ids_turquoise) #232
intersect(deg_list$id, ids_brown) #5
intersect(deg_list$id,ids_grey) #155

comumn_genes =intersect(deg_list$id, colnames(datExpr))


############################
TOM = TOMsimilarityFromExpr(datExpr, power = 9);
# Select modules
module_t = c("turquoise");
# Select module probes
probes = colnames(datExpr)
inModule = is.finite(match(moduleColors, module_t));
modProbes = probes[inModule];
modGenes = annotation$external_gene_name[match(modProbes, annotation$ensembl_transcript_id)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
modTOM <- as.data.frame(modTOM)

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

chooseTopHubInEachModule(datExpr,moduleColors,power = 9,type = "signed")

n_occur <- data.frame(table(nodes$Label))
n_occur[n_occur$Freq > 1,]

write.csv(edges,file = "160_wgcna/edges.csv", row.names = F)
write.csv(nodes,file = "160_wgcna/nodes.csv", row.names = F)



##############GO E KEGG################
# Select modules
modules = c("turquoise");
valid_entrez = annotation$entrezgene_id[match(valid_ids, annotation$ensembl_transcript_id)]
valid_ids
GO_results <- enrichGO(valid_entrez, "org.Hs.eg.db", ont = "BP")


fit = plot(barplot(GO_results, showCategory = 10, main = "GO"))
png("out_turquoise.png", res = 250,width = 1400, height = 2000)
print(fit)
dev.off()

###########################################
###################KEGGREST####################
valid_symbols_turquoise = annotation$entrezgene_id[match(valid_ids, annotation$ensembl_transcript_id)]
keg_turquoise = enrichKEGG(valid_symbols_turquoise, keyType = "kegg", pvalueCutoff=0.5)
print(keg_turquoise)
kegg_t = plot(barplot(keg_turquoise, showCategory = 10))
png("kegg_turquoise.png", res = 250, width = 1400, height = 2000)
print(kegg_t)
dev.off()

# Recalculate topological overlap if needed
############################
TOM = TOMsimilarityFromExpr(datExpr, power = 9);
# Select modules
module_blue = c("blue");
# Select module probes
probes = colnames(datExpr)
inModule_blue = is.finite(match(moduleColors, module_blue));
modProbes_blue = probes[inModule_blue];
modGenes_blue = annotation$external_gene_name[match(modProbes_blue, annotation$ensembl_transcript_id)];
# Select the corresponding Topological Overlap
modTOM_blue = TOM[inModule_blue, inModule_blue];
dimnames(modTOM_blue) = list(modProbes_blue, modProbes_blue)
modTOM_blue <- as.data.frame(modTOM_blue)

valid_ids_blue = deg_list$id[deg_list$id %in% colnames(modTOM_blue)]

modTOM_filter_blue = modTOM_blue[valid_ids_blue, valid_ids_blue]

cyt_blue = exportNetworkToCytoscape(modTOM_filter_blue,
                               edgeFile = paste("edges-", paste(module_blue, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("nodes-", paste(module_blue, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.1, #0,05 0,1 0.15
                               nodeNames = modProbes_blue,
                               altNodeNames = modGenes_blue,
                               nodeAttr = moduleColors[inModule_blue])
edges_blue <- cyt_blue$edgeData
colnames(edges_blue)[1] <- "Source"
colnames(edges_blue)[2] <- "Target"
edges_blue <-edges_blue[, -4:-6]

nodes_blue <- cyt_blue$nodeData
colnames(nodes_blue)[1] <- "id"
colnames(nodes_blue)[2] <- "Label"

chooseTopHubInEachModule(datExpr,moduleColors,power = 9,type = "signed")

n_occur_blue <- data.frame(table(nodes_blue$Label))
n_occur_blue[n_occur_blue$Freq > 1,]

write.csv(edges_blue,file = "160_wgcna/edges_blue.csv", row.names = F)
write.csv(nodes_blue,file = "160_wgcna/nodes_blue.csv", row.names = F)


##############GO E KEGG################
# Select modules
module_blue = c("blue");

valid_entrez_blue = annotation$entrezgene_id[match(valid_ids_blue, annotation$ensembl_transcript_id)]
valid_ids_blue
GO_results_blue <- enrichGO(valid_entrez_blue, "org.Hs.eg.db", ont = "BP")


fit_blue = plot(barplot(GO_results_blue, showCategory = 10, main = "GO"))
png("out_blue.png", res = 250, width = 1400, height = 2000)
print(fit_blue)
dev.off()

valid_symbols_blue = annotation$entrezgene_id[match(valid_ids_blue, annotation$ensembl_transcript_id)]
keg_blue = enrichKEGG(valid_symbols_blue, keyType = "kegg", pvalueCutoff=0.5)
print(keg_blue)
kegg_blue = plot(barplot(keg_blue, showCategory = 10))
png("kegg_blue.png", res = 250, width = 1400, height = 2000)
print(kegg_blue)
dev.off()



