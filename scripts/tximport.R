library(tximport)
library(readr)

# Ler samples
samples <- read.table("samples.tsv", header=TRUE, sep="\t")

# Caminhos dos arquivos quant.sf
files <- file.path("results/salmon", samples$sample_id, "quant.sf")
names(files) <- samples$sample_id

# Ler tx2gene
tx2gene <- read.csv("reference/tx2gene.csv")

# tximport
txi <- tximport(files,
                type = "salmon",
                tx2gene = tx2gene)

# Salvar resultados
write.csv(as.data.frame(txi$counts), "results/deseq2/counts.csv")
write.csv(as.data.frame(txi$length), "results/deseq2/lengths.csv")