#!/usr/bin/env Rscript

suppressMessages(library("DRIMSeq"))
suppressMessages(library("GenomicFeatures"))

cat("Loading counts, conditions and parameters.\n")
cts <- as.matrix(read.csv("merged/all_counts.tsv", sep="\t", row.names="Reference", stringsAsFactors=FALSE))
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

design <- model.matrix(~condition, data=DRIMSeq::samples(d))

suppressMessages(library("dplyr"))

# Sum transcript counts into gene counts:

trs_cts <- counts(d)
gene_cts <- trs_cts %>% select(c(1, 3:ncol(trs_cts)))  %>% group_by(gene_id) %>% summarise_all(funs(sum)) %>% data.frame()
rownames(gene_cts) <- gene_cts$gene_id
gene_cts$gene_id <- NULL

head(trs_cts)
head(gene_cts)

# Differential gene expression using edgeR:
suppressMessages(library("edgeR"))

cat("foo\n")

y <- DGEList(gene_cts)
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
edger_res <- topTags(qlf, n=nrow(y), sort.by="PValue")[[1]]

head(edger_res)

