#!/usr/bin/env Rscript

cat("Loading DESeq2.\n")
suppressMessages(library("DESeq2"))

cat("Loading counts and conditions.\n")
cts <- as.matrix(read.csv("merged/all_counts.tsv", sep="\t", row.names="Reference"))
coldata <- read.csv("de_analysis/coldata.tsv", row.names="sample", sep="\t")

cat("Performing differential expression analysis.\n")
cts <- cts[, rownames(coldata)]
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "untreated")
dds <- DESeq(dds)
res <- results(dds, alpha=0.05)
resOrdered <- res[order(res$pvalue),]

summary(resOrdered)
write.table(as.data.frame(resOrdered), file="de_analysis/deseq2_results.tsv", sep="\t")

cat("Generating plots.\n")
suppressMessages(library("ggplot2"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("vsn"))
suppressMessages(library("pheatmap"))
suppressMessages(library("hexbin"))

pdf("de_analysis/deseq2_plots.pdf")
plotMA(resOrdered, ylim=c(-3,3), main="MA plot")

ggplot(data.frame(resOrdered), aes(x=log2FoldChange, 
                        y=-log10(resOrdered$pvalue))) + 
                        geom_point(alpha = 0.6) +
                        theme_bw() +
                        xlab("Log 2 Fold Change") +
                        ylab("-log10 Adjusted pvalue") +
                        ggtitle("Volcano plot: Control vs. Treated")


vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, colnames(vsd), sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(vsd, intgroup=c("condition"))

plotDispEsts(dds)

meanSdPlot(assay(vsd))

