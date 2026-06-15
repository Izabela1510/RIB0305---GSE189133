import pandas as pd

samples = pd.read_csv("samples.tsv", sep="\t")

# Expandir SRRs por amostra
srr_list = []
for _, row in samples.iterrows():
    for srr in row["srr"].split(","):
        srr_list.append(srr.strip())

rule all:
    input:
        expand("results/salmon/{srr}/quant.sf", srr=srr_list),
        "results/multiqc/multiqc_report.html",
        "results/deseq2/contrast_results.csv",
        "results/enrichment/go_enrichment.csv"

rule salmon_quant:
    input:
        fq = "data/raw/{srr}.fastq.gz",
        index = "reference/salmon_index"
    output:
        "results/salmon/{srr}/quant.sf"
    threads: 4
    shell:
        """
        salmon quant -i {input.index} -l A -r {input.fq} \
            -p {threads} --gcBias -o results/salmon/{wildcards.srr}
        """

rule tximport:
    input:
        expand("results/salmon/{srr}/quant.sf", srr=srr_list),
        tx2gene = "reference/tx2gene.csv",
        samples = "samples.tsv"
    output:
        counts = "results/deseq2/counts.csv"
    script:
        "scripts/tximport.R"

rule deseq2:
    input:
        counts = "results/deseq2/counts.csv",
        samples = "samples.tsv"
    output:
        "results/deseq2/contrast_results.csv"
    script:
        "scripts/deseq2.R"

rule enrichment:
    input:
        "results/deseq2/contrast_results.csv"
    output:
        "results/enrichment/go_enrichment.csv"
    script:
        "scripts/enrichment.R"
