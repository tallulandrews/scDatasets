# Tusi et al. Population snapshots predict early haematopoietic and erythroid hierarchies. Nature. 2018. 
#
# PMID: 29466336
# GSE89754
files_raw <- c("GSM2388072_basal_bone_marrow.raw_umifm_counts.csv.gz", "GSM2388073_epo_bone_marrow.raw_umifm_counts.csv.gz", "GSM2388074_fetal_liver.raw_umifm_counts.csv.gz")
files_P_raw <- c("GSM2985844_P1.raw_umifm_counts.csv.gz", "GSM2985845_P1-CD71hi.raw_umifm_counts.csv.gz", "GSM2985846_P2.raw_umifm_counts.csv.gz", "GSM2985847_P3.raw_umifm_counts.csv.gz", "GSM2985848_P4.raw_umifm_counts.csv.gz", "GSM2985849_P5.raw_umifm_counts.csv.gz")
files_norm <- c("GSM2388072_basal_bone_marrow.filtered_normalized_counts.csv.gz", "GSM2388073_epo_bone_marrow.filtered_normalized_counts.csv.gz", "GSM2388074_fetal_liver.filtered_normalized_counts.csv.gz")


count_mat <- vector()
anno <- vector();

for (f in files_raw) {
	print(f)
	require("Matrix")
	mat <- read.table(f, sep=",", header=T)
	g_i <- min(which(grepl("Rik", colnames(mat))))-1
	a <- mat[,1:g_i]
	mat <- mat[,-c(1:g_i)]
	mat <- t(mat)
	genes <- rownames(mat);
	genes <- sub("^X([1234567890]+)", "\\1", genes);
	rownames(mat) <- genes
	sample_names <- a[,1]
	if (length(unique(sample_names)) < nrow(a)) {
		sample_names <- paste(as.character(sample_names), 1:nrow(a), sep="")
	}
	colnames(mat) <- sample_names
	rownames(a) <- sample_names
	mat <- Matrix(mat)
	count_mat <- cbind(count_mat, mat)
	anno <- rbind(anno, a)
}

require("SingleCellExperiment")
obj <- SingleCellExperiment(assays=list(counts=count_mat), colData=anno)
saveRDS(obj, "Tusi_rawCounts.rds")

obj_clean <- obj[,obj$pass_filter==1]

require("scater")
obj_clean <- normalize(obj_clean)
obj <- obj_clean

markers <- read.table("/nfs/users/nfs_t/ta6/Data/Hematopoiesis_Markers.txt", header=TRUE)
markers <- markers[order(markers[,1]),]

obj_clean <- obj_clean[Matrix::rowMeans(obj_clean@assays[["counts"]]>0) > 0.01,]

# 50 PCs - irlba
require(irlba)
pcs <- irlba(obj_clean@assays[["logcounts"]], nv=50)

# k = 40 - Rphenograph
require(Rphenograph)
clusters <- Rphenograph::Rphenograph(pcs$v, 40)

#obj_clean <- obj_clean[rownames(obj_clean) %in% markers[,1],]
#clusters <- kmeans(as.matrix(t(obj_clean@assays[["logcounts"]])), centers=10)

row_mean_aggregate <- function (mat, groups) {
    MAT <- as.matrix(mat)
    x <- split(seq(ncol(MAT)), groups)
    result <- sapply(x, function(a) Matrix::rowMeans(MAT[, a]))
    return(result)
}

profiles <- row_mean_aggregate(obj_clean@assays[["logcounts"]], clusters[[2]]$membership)
#profiles <- row_mean_aggregate(obj_clean@assays[["logcounts"]], clusters$cluster)
rownames(profiles) <- rownames(obj_clean)
profiles <- t(apply(profiles, 1, scale))

tmp <- vector();
for (i in as.character(unique(markers[,2]))) {
        score<- colMeans(profiles[rownames(profiles) %in% as.character(markers[markers[,2]==i,1]),])
	tmp <- rbind(tmp, score)
}
thresh <- min(round(apply(tmp,1,max)))
tmp[tmp < thresh] <- 0

groups <- as.character( clusters[[2]]$membership )
for (i in 1:ncol(tmp)) {
	if (sum(tmp[,i]) ==0) {next;}
	g <- which(tmp[,i] == max(tmp[,i]))
	groups[groups == i] <- rownames(tmp)[g]
}

#groups[groups %in% as.character( clusters[[2]]$membership )] <- paste("Undiff", groups[groups %in% as.character( clusters[[2]]$membership )], sep="")

obj$type <- groups;
saveRDS(obj, "Tusi_clean.rds")
