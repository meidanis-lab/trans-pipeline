BiocManager::install("GWENA")
library("GWENA")
library("biomaRt")
threads_to_use <- 2

#Reading counts
counts = read.csv("100_counts/counts.csv", header = T, row.names = 1, sep = ",")
meta = read.csv("105_meta/meta.csv", header=T, row.names=1)


######################################Normalizing
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~Status)

vsd <- vst(dds, blind = T)
datExpr0 = t(assay(vsd))
datExpr = filter_RNA_seq(datExpr0, min_count = 5, method = "all")
min(datExpr)
datExpr = t(datExpr)

# Calculate the variance or standard deviation for each gene
gene_variability <- apply(datExpr, 1, var)  # Using variance

# Determine the number of genes to select (5000 most variant genes)
num_genes <- round(0.50 * length(gene_variability))

# Sort the genes based on their variability in descending order
sorted_genes <- order(gene_variability, decreasing = TRUE)

# Select the top 5000 genes with the highest variability
selected_genes <- sorted_genes[1:num_genes]

# Subset the datExpr matrix using the selected genes
datExpr <- datExpr[selected_genes,]
datExpr_filtered = t(datExpr)


#Building net
net <- build_net(datExpr_filtered, cor_func = "spearman", 
                 n_threads = 5)
net$metadata$power
fit_power_table <- net$metadata$fit_power_table
fit_power_table[fit_power_table$Power == net$metadata$power, "SFT.R.sq"]


##Detecting modules
modules <- detect_modules(datExpr_filtered, 
                          net$network, 
                          detailled_result = TRUE,
                          merge_threshold = 0.25)
# Number of modules before merging :
length(unique(modules$modules_premerge))
# Number of modules after merging: 
length(unique(modules$modules))

#Visualization of modules
layout_mod_merge <- plot_modules_merge(
  modules_premerge = modules$modules_premerge, 
  modules_merged = modules$modules)

ggplot2::ggplot(data.frame(modules$modules %>% stack), 
                ggplot2::aes(x = ind)) + ggplot2::stat_count() +
  ggplot2::ylab("Number of genes") +
  ggplot2::xlab("Module")

#RElating modules with functional enchriment
enrichment <- bio_enrich(modules$modules)
plot_enrichment(enrichment)

#Association with phenotype
mapping <- list(Metastasis = c("no" = 0, "Lung" = 1, "Lung_bones" = 2),
                Status = c("Dead" = 1, "Alive" = 0),
                TP53 = c("no" = 0, "Mut" = 1),
                COPS3 = c("no" = 0, "COP" = 1),
                RB1 = c("no" = 0,"Mut" = 1),
                CDK = c("no" = 0, "Mut" = 1),
                PI3K.mTOR = c("no" = 0, "PTEN_Mut" = 1, "PTEN_del" = 2, 
                              "NF1_Mut" = 3, "NF1_del" = 4, "AKT1_Mut" = 5,
                              "NF1_Mut.PDPK1_Mut" = 6))

# Define a function to convert categorical values to numerical values
for (col in colnames(meta)) {
  if (col %in% names(mapping)) {
    meta[[col]] <- mapping[[col]][meta[[col]]]
  }
}

phenotype_association <- associate_phenotype(
  modules$modules_eigengenes, 
  meta %>% dplyr::select(Metastasis,Status,TP53,RB1,CDK,COPS3,PI3K.mTOR))

plot_modules_phenotype(phenotype_association)

association = phenotype_association[["association"]][phenotype_association[["association"]] >= 0.1 | phenotype_association[["association"]] <= -0.1, ]
association = association[complete.cases(association),]



#####Annotation with Biomart
ensembl_id = as.data.frame(colnames(datExpr_filtered))

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
                   values = ensembl_id$`colnames(datExpr_filtered)`,
                   mart = ensembl_con)
library(igraph)
module1 <- modules$modules$`1`
mod1 = intersect(module1, deg_list$id)
enrichment1 <- bio_enrich(mod1)
plot_enrichment(enrichment1)
probes = annotation$external_gene_name[match(mod1, annotation$ensembl_transcript_id)]
adj_matrix_1 = net$network[mod1, mod1]
dimnames(adj_matrix_1) = list(probes, probes)
write.csv(adj_matrix_1, file = "205_gwena/adj_1.csv")


module2 = modules$modules$'2'
mod2 =intersect(module2, deg_list$id)
enrichment2 <- bio_enrich(mod2)
plot_enrichment(enrichment2)
probes = annotation$external_gene_name[match(mod2, annotation$ensembl_transcript_id)]
adj_matrix_2 = net$network[mod2, mod2]
dimnames(adj_matrix_2) = list(probes, probes)
write.csv(adj_matrix_2, file = "205_gwena/adj_2.csv")


module3 = modules$modules$`3`
mod3 =intersect(module3, deg_list$id)
enrichment3 <- bio_enrich(mod3)
plot_enrichment(enrichment3)
probes = annotation$external_gene_name[match(mod3, annotation$ensembl_transcript_id)]
adj_matrix_3 = net$network[mod3, mod3]
dimnames(adj_matrix_3) = list(probes, probes)
write.csv(adj_matrix_3, file = "205_gwena/adj_3.csv")


module4 = modules$modules$`4`
mod4 =intersect(module4, deg_list$id)
enrichment4 <- bio_enrich(mod4)
plot_enrichment(enrichment4)
probes = annotation$external_gene_name[match(mod4, annotation$ensembl_transcript_id)]
adj_matrix_4 = net$network[mod4, mod4]
dimnames(adj_matrix_4) = list(probes, probes)
write.csv(adj_matrix_4, file = "205_gwena/adj_4.csv")


