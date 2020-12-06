data = read.table("GSE72857_umitab.txt.gz")
design = read.delim("GSE72857_experimental_design.txt")
#anno = read.delim("GSE72857_series_matrix.txt.gz")

source("~/R-Scripts/Ensembl_Stuff.R")
musg_genes <- musg_name_map[,2]

source("~/Data_Processing_Scripts/Fix_Gene_Names.R")
rownames(data) <- fix_gene_names(rownames(data), musg_genes)

rownames(design) <- design[,1]
design <- design[match(colnames(data), rownames(design)),]
design$facs_type <- rep("None", nrow(design))
design$facs_type[design$CD34_protein < 200 & design$FcgR3_protein < 100] <- "MEP"
design$facs_type[design$CD34_protein > 200 & design$FcgR3_protein > 100 & design$FcgR3_protein < 600] <- "CMP"
design$facs_type[design$CD34_protein > 200 & design$FcgR3_protein > 600] <- "GMP"

require("SingleCellExperiment")
require("Matrix")
require("scater")
obj <- SingleCellExperiment(assays=list(counts=Matrix(as.matrix(data))), colData=design)
obj <- normalize(obj)
saveRDS(obj, "Paul.rds")

obj <- obj[,Matrix::colSums(obj@assays[["logcounts"]]>0) >1000]
obj <- obj[,obj$Batch_desc %in% c("Unsorted myeloid", "CMP CD41")]

set.seed(2010)
obj <- plotPCA(obj, return_SCE=TRUE)

marker_set <- c("Ccl5", "Prg2", "Prss34", "Meis1", "H2-Aa", "Cd74", "Pf4", "Pbx1", "Serpina3f", "Apoe", "Gata2", "Gata1", "Gfi1b", "Car1", "Car2", "Klf1", "Zfpm1", "Cpox", "Beta-s", "Hbb-b1", "Hba-a2", "Cebpe", "Csf1r", "Pu.1", "Lgals1", "Irf8", "Elane", "Prtn3", "Mpo", "Flt3", "Ifitm1", "Lmo4")

set.seed(2010)
thing <- kmeans(t(obj@assays[["logcounts"]][rownames(obj) %in% marker_set,]), centers=4, nstart=50)
colData(obj)$cluster <- thing$cluster
obj <- plotPCA(obj, return_SCE=TRUE, feature_set=marker_set, colour_by="cluster", shape_by="facs_type")
a <- table(obj$cluster, obj$Batch_desc)
b <- table(obj$cluster, obj$facs_type)
cmp_cluster <- which(b[,1] == max(b[,1]))
gmp_cluster <- which(b[,2] == max(b[,2]))
mep_cluster <- which(b[,3] == max(b[,3]))
cmpCD41_cluster <- which(a[,5] == max(a[,5]))

colData(obj)$named_clusters <- rep("outliers", times=ncol(obj))

obj$named_clusters[obj$cluster == cmp_cluster] <- "CMP"
obj$named_clusters[obj$cluster == cmpCD41_cluster] <- "CMP-CD41"
obj$named_clusters[obj$cluster == mep_cluster] <- "MEP"
obj$named_clusters[obj$cluster == gmp_cluster] <- "GMP"
saveRDS(obj, "Paul_clean.rds")
