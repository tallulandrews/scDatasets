# Mohammed et al. 2017 Single-Cell Landscape of Transcriptional Heterogeneity and Cell Fate Decisions during Mouse Early Gastrulation. Cell Reports. 
# PMID: 28768204
# GSE100597

a <- read.table("GSE100597_count_table_QC_filtered.txt.gz")
ids <- strsplit(colnames(a), "_")
time <- unlist(lapply(ids, function(x)x[1]))
num <- unlist(lapply(ids, function(x)x[length(x)]))
sample <- unlist(lapply(ids, function(x)x[length(x)-1]))

require("scater")
require("Matrix")
obj <- SingleCellExperiment(assays=list(counts=Matrix(as.matrix(a))), colData=data.frame(Time=time, numCells=num, Sample=sample, cell_type1=time))
saveRDS(obj, file="Gastrulation2.rds")
