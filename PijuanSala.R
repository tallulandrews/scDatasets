# Pijuan-Sala et al. A single-cell molecular map of mouse gastrulation and early organogenesis. Nature. 2019
# PMID: 
# https://github.com/MarioniLab/EmbryoTimecourse2018

# 10X
# QC - >=1,000 expressed genes, mito < 2%, mean log2(norm) >= 10^-3, exclude chrY, Xist, tdTomato construct
# Normalization: scran - quickCluster(method=igraph, min=100, max=3000)
# HVG: scran - trendVar and decomposeVar, loess span 0.05, FDR 5%
# Batch effects: scran - fastMNN on 50 PCs from HVGs
# doublet removal: scan - doubletCells, 50 PCs, 10-NN (buildSNNGraph in scran), louvain clustering from igraph, hierarchically applied, FDR 10%

# SS2
# nuclear reads > 50,000, genes > 4000, mito < 10%

# Clustering: buildSSNGraph (scran) 50 batch-corrected PCs on HVGs, 10-NN
#		cluster_louvain (igraph) x2 levels 

require("Matrix")
require("SingleCellExperiment")
anno <- read.table("atlas/meta.tab", sep="\t", header=T)
count_mat <- readMM("atlas/raw_counts.mtx")
count_mat <- as(count_mat, "dgCMatrix")
gene_ids <- read.table("atlas/genes.tsv", stringsAsFactors=FALSE)
cell_ids <- read.table("atlas/barcodes.tsv", stringsAsFactors=FALSE)
pcs <- readRDS("atlas/corrected_pcas.rds")
sfs <- read.table("atlas/sizefactors.tab")
rownames(count_mat) <- gene_ids[,1]
colnames(count_mat) <- cell_ids[,1]
rownames(anno) <- anno[,1]
rownames(gene_ids) <- gene_ids[,1]
colnames(gene_ids) <- c("ensg", "feature_symbol")

sce <- SingleCellExperiment(assays=list(counts=count_mat), colData=anno, rowData=gene_ids)
sce@reducedDims <- SimpleList(pca=pcs$all)
sce@int_colData <- DataFrame(size_factor=sfs[,1])

require("scater")
sce <- normalize(sce)
saveRDS(sce, "PijuanSala.rds")


require("Matrix")
require("SingleCellExperiment")
anno <- read.table("atlas/meta.tab", sep="\t", header=T)
keep <- !is.na(anno$celltype)
anno <- anno[keep,]
count_mat <- readMM("atlas/raw_counts.mtx")
count_mat <- as(count_mat[,keep], "dgCMatrix")
gene_ids <- read.table("atlas/genes.tsv", stringsAsFactors=FALSE)
cell_ids <- read.table("atlas/barcodes.tsv", stringsAsFactors=FALSE)
pcs <- readRDS("atlas/corrected_pcas.rds")
sfs <- read.table("atlas/sizefactors.tab")
sfs <- sfs[,keep]
pcs <- pcs[,keep]
rownames(count_mat) <- gene_ids[,1]
colnames(count_mat) <- cell_ids[keep,1]
rownames(anno) <- anno[,1]
rownames(gene_ids) <- gene_ids[,1]
colnames(gene_ids) <- c("ensg", "feature_symbol")

sce <- SingleCellExperiment(assays=list(counts=count_mat), colData=anno, rowData=gene_ids)
sce@reducedDims <- SimpleList(pca=pcs$all)
sce@int_colData <- DataFrame(size_factor=sfs[,1])

require("scater")
sce <- normalize(sce)
saveRDS(sce, "PijuanSala_clean.rds")


