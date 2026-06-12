library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# Ler resultados DESeq2
res <- read.csv("results/deseq2/contrast_results.csv", row.names=1)

# Filtrar genes DE
de_genes <- rownames(res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, ])

# Converter para Entrez ID
entrez <- bitr(de_genes, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

# Enriquecimento GO
go_enrich <- enrichGO(
  gene = entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)
write.csv(as.data.frame(go_enrich), "results/enrichment/go_enrichment.csv")
dotplot(go_enrich, showCategory=20)
ggsave("results/enrichment/go_plot.png")

# Enriquecimento KEGG
kegg_enrich <- enrichKEGG(
  gene = entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)
write.csv(as.data.frame(kegg_enrich), "results/enrichment/kegg_enrichment.csv")
dotplot(kegg_enrich, showCategory=20)
ggsave("results/enrichment/kegg_plot.png")