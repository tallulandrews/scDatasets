#Rodriguez-Fraticelli et al. Clonal analysis of lineage fate in native haematopoiesis. Nature 2018
# PMID: 29323290
# GEO: GSE90742
set.seed(1910)

files <- c("GSM2411664_LTHSC.raw_umifm_counts.csv.gz", "GSM2411665_MPP2.raw_umifm_counts.csv.gz", "GSM2411666_MPP3.raw_umifm_counts.csv.gz", "GSM2411667_MPP4.raw_umifm_counts.csv.gz", "GSM2411668_STHSC.raw_umifm_counts.csv.gz")
require("Matrix")
expr_mat <- vector()
anno <- vector()
for (f in files) {
	dat <- read.table(f, sep=",", header=T, stringsAsFactors=FALSE)
	rownames(dat) <- dat[,1]
	ann <- dat[,1:5]
	dat <- dat[,-c(1:5)]
	genes <- colnames(dat);
	dat <- t(dat);
	genes <- sub("^X([1234567890]+)","\\1", genes)
	rownames(dat) <- genes;
	dat <- Matrix(dat)
	expr_mat <- cbind(expr_mat, dat)
	anno <- rbind(anno, ann)
}

require("SingleCellExperiment")
obj <- SingleCellExperiment(assays=list(counts=expr_mat), colData=anno);

saveRDS(obj, "Rodriguez_raw.rds")

require("scater")
obj <- obj[,colData(obj)$pass_filter == 1]
lib.size <- Matrix::colSums(obj@assays[["counts"]])
obj@assays[["norm"]] <- t( t(obj@assays[["counts"]])/lib.size*median(lib.size) )
# gene filter cv > 2, & mean > 0.05
g_means <- Matrix::rowMeans(obj@assays[["norm"]])
g_var <- Matrix::rowSums(( obj@assays[["norm"]]-g_means )^2) /(ncol(obj)-1)
g_cv <- sqrt(g_var)/g_means
g_cv[is.na(g_cv)] <- 0
gene_filter <- g_means > 0.05 & g_cv > 2
obj_clean <- obj[gene_filter,]

require("CycleMix")
cc_genes <- as.character(MGeneSets[[1]][,1])
obj_clean <- obj_clean[! rownames(obj_clean) %in% cc_genes,]

# 50 PCs - irlba
require(irlba)
obj_clean <- scater::normalize(obj_clean);
pcs <- irlba(obj_clean@assays[["logcounts"]], nv=50)

# k = 40 - Rphenograph
require(Rphenograph)
clusters <- Rphenograph::Rphenograph(pcs$v, 40)


# match to names
marks <- read.delim("Rodriguez_markers.csv", header=T, sep="\t")
row_mean_aggregate <- function (mat, groups) {
    MAT <- as.matrix(mat)
    x <- split(seq(ncol(MAT)), groups)
    result <- sapply(x, function(a) Matrix::rowMeans(MAT[, a]))
    return(result)
}

profiles <- row_mean_aggregate(obj_clean@assays[["logcounts"]], clusters[[2]]$membership)
rownames(profiles) <- rownames(obj_clean)
groups <- as.character( clusters[[2]]$membership )
for (i in as.character(unique(marks[,6]))) {
	score<- colMeans(profiles[rownames(profiles) %in% as.character(marks[marks[,6]==i,1]),])
	id <- which(score == max(score))
	groups[groups == as.character(id)] <- i
}
groups[groups %in% as.character( clusters[[2]]$membership )] <- paste("Undiff", groups[groups %in% as.character( clusters[[2]]$membership )], sep="")

names(groups) <- colnames(obj_clean)
obj$cluster <- groups
saveRDS(obj, "Rodriguez_clean.rds")
