library(tximport)
library(readr)

# Ler samples
samples <- read.table("samples.tsv", header=TRUE, sep="\t")

# Expandir: uma linha por SRR, mantendo referência à amostra
srr_to_sample <- data.frame()
for (i in seq_len(nrow(samples))) {
  srrs <- trimws(strsplit(samples$srr[i], ",")[[1]])
  for (srr in srrs) {
    srr_to_sample <- rbind(srr_to_sample, data.frame(srr=srr, sample_id=samples$sample_id[i]))
  }
}

# Arquivos quant.sf por SRR
files <- file.path("results/salmon", srr_to_sample$srr, "quant.sf")
names(files) <- srr_to_sample$srr

# Ler tx2gene
tx2gene <- read.csv("reference/tx2gene.csv")

# tximport por SRR
txi <- tximport(files,
                type = "salmon",
                tx2gene = tx2gene,
                ignoreTxVersion = TRUE)

# Agregar réplicas técnicas (somar counts por amostra)
sample_ids <- unique(srr_to_sample$sample_id)
counts_agg <- sapply(sample_ids, function(s) {
  srrs <- srr_to_sample$srr[srr_to_sample$sample_id == s]
  rowSums(txi$counts[, srrs, drop=FALSE])
})

# Criar diretório
dir.create("results/deseq2", recursive=TRUE, showWarnings=FALSE)

# Salvar
saveRDS(txi, "results/deseq2/txi.rds")
write.csv(as.data.frame(counts_agg), "results/deseq2/counts.csv")

cat("tximport concluído!\n")
