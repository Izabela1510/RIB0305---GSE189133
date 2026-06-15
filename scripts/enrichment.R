library(topGO)
library(org.Hs.eg.db)
library(ggplot2)

# Ler resultados DESeq2
res <- read.csv("results/deseq2/contrast_results.csv", row.names=1)

# Genes DE
de_genes <- rownames(res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ])

# Gene list para topGO (1 = DE, 0 = não DE)
all_genes <- rownames(res)
gene_list <- factor(as.integer(all_genes %in% de_genes))
names(gene_list) <- all_genes

# topGO object
go_data <- new("topGOdata",
  ontology = "BP",
  allGenes = gene_list,
  annot = annFUN.org,
  mapping = "org.Hs.eg.db",
  ID = "ensembl")

# Teste
result_fisher <- runTest(go_data, algorithm="classic", statistic="fisher")
result_elim <- runTest(go_data, algorithm="elim", statistic="fisher")

# Tabela de resultados
go_table <- GenTable(go_data,
  classicFisher = result_fisher,
  elimFisher = result_elim,
  orderBy = "classicFisher",
  topNodes = 50)

dir.create("results/enrichment", recursive=TRUE, showWarnings=FALSE)
write.csv(go_table, "results/enrichment/go_enrichment.csv", row.names=FALSE)

# Plot top 20 GO terms
go_plot <- go_table[1:20, ]
go_plot$pval <- as.numeric(go_plot$classicFisher)
go_plot$Term <- factor(go_plot$Term, levels=rev(go_plot$Term))

ggplot(go_plot, aes(x=-log10(pval), y=Term, size=Significant)) +
  geom_point(color="steelblue") +
  theme_bw() +
  labs(title="Top 20 GO Biological Process", x="-log10(p-value)", y="GO Term")

ggsave("results/enrichment/go_plot.png", width=10, height=8)

cat("Enriquecimento concluído!\n")
cat("GO terms significativos:", sum(as.numeric(go_table$classicFisher) < 0.05, na.rm=TRUE), "\n")
