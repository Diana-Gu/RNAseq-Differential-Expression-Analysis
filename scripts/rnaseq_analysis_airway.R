# RNA-seq differential expression analysis using real dataset
# Dataset: airway (RNA-seq data of dexamethasone treatment)

install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("airway")
BiocManager::install("pheatmap")
BiocManager::install("matrixStats")

# Load libraries
library(DESeq2)
library(airway)
library(ggplot2)
library(pheatmap)
library(matrixStats)

# Load airway dataset
data(airway)

# Convert to DESeq2 dataset
dds <- DESeqDataSet(airway, design = ~ dex)

# Filter low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run differential expression
dds <- DESeq(dds)

# Extract results
res <- results(dds)

# Save results
write.csv(as.data.frame(res), "results/deseq2_results.csv")

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# PCA plot
pcaData <- plotPCA(vsd, intgroup = "dex", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(PC1, PC2, color = dex)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of RNA-seq samples")

ggsave("results/pca_plot.png", plot = p)

# Heatmap of top variable genes
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 30)
mat <- assay(vsd)[topVarGenes, ]
mat <- mat - rowMeans(mat)

png("results/heatmap_top_genes.png", width = 800, height = 700)
pheatmap(mat, annotation_col = as.data.frame(colData(dds)))
dev.off()

# Volcano plot
res_df <- as.data.frame(res)

res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "yes", "no")

v <- ggplot(res_df, aes(log2FoldChange, -log10(padj), color = significant)) +
  geom_point(alpha = 0.7) +
  ggtitle("Volcano plot of differential expression")

ggsave("results/volcano_plot.png", plot = v)
