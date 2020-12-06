# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77288 
# Tung et al. (2017) Batch effects and the effective design of single-cell gene expression studies. Sci Rep. (7):39921
# PMID: 28045081
# GEO: GSE77288

reads <- read.table("GSE77288_reads-raw-single-per-sample.txt.gz", header=T)
umi <- read.table("GSE77288_molecules-raw-single-per-sample.txt.gz", header=T)

#anno <- read.table("GSE77288.txt");
anno <- umi[,1:3]
umi <- umi[,-c(1:3)]
umi <- as.matrix(t(umi))
reads <- as.matrix(t(reads[,-c(1:3)]));

cellIDs <- paste("cell", 1:ncol(umi), sep="")
colnames(umi) <- cellIDs
colnames(reads) <- cellIDs
rownames(anno) <- cellIDs

anno$cell_type1 <- anno$individual

source("~/R-Scripts/Ensembl_Stuff.R")
gene_info <- data.frame(ensg=rownames(umi), feature_symbol=General_Map(rownames(umi), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol"))
rownames(gene_info) <- rownames(umi)

require("SingleCellExperiment")

sce <- SingleCellExperiment(assays=list(counts=umi, reads=reads), colData=anno, rowData=gene_info);

saveRDS(sce, "tung.rds")
