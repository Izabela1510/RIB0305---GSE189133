library(DESeq2)
library(ggplot2)
library(pheatmap)
library(readr)

# Ler dados
counts <- read.csv("results/deseq2/counts.csv", row.names=1)
samples <- read.table("samples.tsv", header=TRUE, sep="\t")

# Criar objeto DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData = samples,
  design = ~ fam193a * nutlin
)

# Filtrar genes com baixa expressão
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Rodar DESeq2
dds <- DESeq(dds)

# Resultados
res <- results(dds, alpha = 0.05)
res_df <- as.data.frame(res)
write.csv(res_df, "results/deseq2/contrast_results.csv")

# PCA plot
vsd <- vst(dds, blind=FALSE)
pca_data <- plotPCA(vsd, intgroup=c("fam193a","nutlin"), returnData=TRUE)
ggplot(pca_data, aes(PC1, PC2, color=fam193a, shape=nutlin)) +
  geom_point(size=3) +
  theme_bw()
ggsave("results/deseq2/pca_plot.png")

# Volcano plot
res_df$sig <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "DE", "NS")
ggplot(res_df, aes(log2FoldChange, -log10(padj), color=sig)) +
  geom_point(alpha=0.5) +
  scale_color_manual(values=c("DE"="red","NS"="grey")) +
  theme_bw()
ggsave("results/deseq2/volcano_plot.png")

# MA plot
png("results/deseq2/ma_plot.png")
plotMA(res)
dev.off()

# Heatmap top 50 genes
top50 <- head(order(res_df$padj), 50)
pheatmap(assay(vsd)[top50,],
         annotation_col=as.data.frame(colData(dds)[,c("fam193a","nutlin")]),
         filename="results/deseq2/heatmap.png")