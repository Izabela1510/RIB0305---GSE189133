# RIB0305---GSE189133
## RNA-seq pipeline: FAM193A KO + Nutlin em células CHP212 de neuroblastoma

### Descrição

Pipeline de análise de RNA-seq para o dataset GSE189133, investigando o efeito do knockout de FAM193A combinado com tratamento por Nutlin-3a em células CHP212 de neuroblastoma humano. Design experimental fatorial 2×2.

## Amostras

| sample_id | fam193a | nutlin | réplica |
|---|---|---|---|
| CHP212_WT_DMSO_1 | WT | control | 1 |
| CHP212_WT_DMSO_2 | WT | control | 2 |
| CHP212_WT_Nutlin_1 | WT | nutlin | 1 |
| CHP212_WT_Nutlin_2 | WT | nutlin | 2 |
| CHP212_KO_DMSO_1 | FAM193AKO | control | 1 |
| CHP212_KO_DMSO_2 | FAM193AKO | control | 2 |
| CHP212_KO_Nutlin_1 | FAM193AKO | nutlin | 1 |
| CHP212_KO_Nutlin_2 | FAM193AKO | nutlin | 2 |

Cada amostra é composta por 2 SRRs (réplicas técnicas agregadas).

## Ferramentas

| Ferramenta | Uso |
|---|---|
| FastQC + MultiQC | Controle de qualidade |
| Salmon | Quantificação de transcritos |
| tximport | Importação e agregação de counts |
| DESeq2 | Expressão diferencial |
| topGO | Enriquecimento GO |

## Resultados principais

- **562 genes diferencialmente expressos** (padj < 0.05, |log2FC| > 1)
  - 310 upregulados
  - 252 downregulados
- Design: `~ fam193a * nutlin`
- Organismo: *Homo sapiens* (GRCh38, Ensembl release 111)

## Como reproduzir

salmon index -t reference/transcriptome.fa.gz -i reference/salmon_index -p 4
Rscript scripts/tximport.R
Rscript scripts/deseq2.R
Rscript scripts/enrichment.R

## Dataset

- **GEO:** [GSE189133](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189133)
- **PMID:** 36897777
