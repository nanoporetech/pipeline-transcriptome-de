#!/usr/bin/env Rscript

cat("Loading DESeq2.\n")
suppressMessages(library("DESeq2"))

cat("Loading counts and conditions.\n")
cts <- as.matrix(read.csv("merged/all_counts.tsv", sep="\t", row.names="Reference"))
coldata <- read.csv("de_analysis/coldata.tsv", row.names="sample", sep="\t")

cat("Performing differential expression analysis.")
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, alpha=0.05)
resOrdered <- res[order(res$pvalue),]

head(resOrdered)
write.csv(as.data.frame(resOrdered), file="de_analysis/deseq2_results.tsv", sep="\t")
