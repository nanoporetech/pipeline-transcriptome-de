#!/usr/bin/env Rscript

suppressMessages(library("DRIMSeq"))
suppressMessages(library("GenomicFeatures"))

cat("Loading counts, conditions and parameters.\n")
cts <- as.matrix(read.csv("merged/all_counts.tsv", sep="\t", row.names="Reference"))
coldata <- read.csv("de_analysis/coldata.tsv", row.names="sample", sep="\t")
coldata$sample_id <- rownames(coldata)

de_params <- read.csv("de_analysis/de_params.tsv", sep="\t", stringsAsFactors=FALSE)

txdb <- makeTxDbFromGFF(de_params$Annotation[[1]])
txdf <- select(txdb, keys(txdb,"GENEID"), "TXNAME", "GENEID")
tab <- table(txdf$GENEID)
txdf$ntx<- tab[match(txdf$GENEID, names(tab))]

strip_version<-function(x) {
    tmp<-data.frame(strsplit(x,".", fixed=TRUE), stringsAsFactors=FALSE)
    tmp<-as.vector(tmp[1,])
    colnames(tmp) <- c()
    rownames(tmp) <- c()
    return(tmp)
}

rownames(cts) <- strip_version(rownames(cts))

cts <- cts[rownames(cts) %in% txdf$TXNAME, ] # FIXME!

txdf <- txdf[match(rownames(cts), txdf$TXNAME), ]
rownames(txdf) <- NULL

counts<-data.frame(gene_id=txdf$GENEID, feature_id=txdf$TXNAME, cts)

d <- dmDSdata(counts=counts, samples=coldata)
#
d <- dmFilter(d, min_samps_gene_expr = de_params$min_samps_gene_expr[[1]], min_samps_feature_expr = de_params$min_samps_feature_expr[[1]],
              min_gene_expr = de_params$min_gene_expr[[1]], min_feature_expr = de_params$min_feature_expr[[1]])

print(d)

# Differential gene expression using edgeR:


