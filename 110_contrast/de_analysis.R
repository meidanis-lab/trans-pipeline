args = commandArgs(trailingOnly=TRUE)

### Usage:
### Rscript --vanilla de_analysis.R <counts.csv> <meta.csv> <column> <control> <condition> <output>

countsfile = args[1]
metafile = args[2]
column = args[3]
print(column)
control = args[4]
print(control)
condition = args[5]
print(condition)
output = args[6]

library("DESeq2")

#Reading the files
counts <- read.csv(countsfile, header = T, row.names = 1, sep = ",")
metadata <- read.csv(metafile, header = T, row.names = 1, sep = ",")

#Checking if row names matches column names 
all(rownames(metadata)==colnames(counts))

#Building the dataset from matrix
dds <- DESeqDataSetFromMatrix(countData = counts,
				colData = metadata,
				design = as.formula(paste("~",column)))

#Differential expression analysis
DE <- DESeq(dds)
RE <- results(DE)

#Contrast
res <- results(DE, contrast = c(as.character(column),
				condition,
				control))


res$sig <- as.factor(ifelse(res$padj < 0.05 & 
				abs(res$log2FoldChange) >2 & 
				res$baseMean > 50,
			ifelse(res$log2FoldChange > 2,'up','down'), 'not'))


res = as.data.frame(res)
res <-  res[order(factor(res$sig, levels = c('up', 'down', 'not'))), ]


summary(res$sig)
write.table(res, file =output)
