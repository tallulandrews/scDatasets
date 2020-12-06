# Peterson et al. Multiplexed quantification of proteins and transcripts in single cells
# PMID: 28854175
# GSE100501

files_RNA <- c(
"GSM2685237_mRNA_1_PBMCs_matrix.txt.gz",
"GSM2685238_mRNA_2_PBMCs_matrix.txt.gz",
"GSM2685239_mRNA_3_PBMCs_matrix.txt.gz",
"GSM2685240_mRNA_4_CD3_enriched_matrix.txt.gz",
"GSM2685241_mRNA_5_CD11b_enriched_matrix.txt.gz",
"GSM2685242_mRNA_6_CD19_enriched_matrix.txt.gz",
"GSM2685248_mRNA_1_Donor_1_no_aCD27_matrix.txt.gz",
"GSM2685249_mRNA_2_Donor_1_no_aCD27_matrix.txt.gz",
"GSM2685250_mRNA_3_Donor_1_with_aCD27_matrix.txt.gz",
"GSM2685251_mRNA_4_Donor_1_with_aCD27_matrix.txt.gz",
"GSM2685252_mRNA_5_Donor_2_no_aCD27_matrix.txt.gz",
"GSM2685253_mRNA_6_Donor_2_with_aCD27_matrix.txt.gz",
"GSM2685254_mRNA_7_Donor_2_with_aCD27_matrix.txt.gz",
"GSM2685255_mRNA_8_Donor_3_no_aCD27_matrix.txt.gz",
"GSM2685256_mRNA_9_Donor_3_with_aCD27_matrix.txt.gz",
"GSM2685257_mRNA_10_Donor_3_with_aCD27_matrix.txt.gz"
)

files_protein <- c(
"GSM2685243_protein_2_PBMCs_matrix.txt.gz",
"GSM2685244_protein_3_PBMCs_matrix.txt.gz",
"GSM2685245_protein_4_CD3_enriched_matrix.txt.gz",
"GSM2685246_protein_5_CD11b_enriched_matrix.txt.gz",
"GSM2685247_protein_6_CD19_enriched_matrix.txt.gz",
"GSM2685258_protein_1_Donor_1_no_aCD27_matrix.txt.gz",
"GSM2685259_protein_2_Donor_1_no_aCD27_matrix.txt.gz",
"GSM2685260_protein_3_Donor_1_with_aCD27_matrix.txt.gz",
"GSM2685261_protein_4_Donor_1_with_aCD27_matrix.txt.gz",
"GSM2685262_protein_5_Donor_2_no_aCD27_matrix.txt.gz",
"GSM2685263_protein_6_Donor_2_with_aCD27_matrix.txt.gz",
"GSM2685264_protein_7_Donor_2_with_aCD27_matrix.txt.gz",
"GSM2685265_protein_8_Donor_3_no_aCD27_matrix.txt.gz",
"GSM2685266_protein_9_Donor_3_with_aCD27_matrix.txt.gz",
"GSM2685267_protein_10_Donor_3_with_aCD27_matrix.txt.gz"
)

count_mat <- vector()
anno <- vector();

for (f in files_RNA) {
	print(f)
	require("Matrix")
	mat <- read.table(f, sep="\t", header=T)
        mat <- Matrix(as.matrix(mat))
	active <- "No"
	if (grepl("with_aCD27", f)) {
		active <- "Yes"
	}
	sample <- sub("_matrix.txt.gz","",f)
	sample <- sub("_with_aCD27","",sample)
	sample <- sub("_no_aCD27","",sample)
	sample <- sub("GSM\\d*_[[:alpha:]]*_","",sample)
	biosample <- sub("^\\d*_","",sample)
	a <- cbind(rep(active, ncol(mat)), rep(sample, ncol(mat)), rep(biosample, ncol(mat)), rep("RNA", ncol(mat)))
	count_mat <- cbind(count_mat, mat)
	anno <- rbind(anno, a)
}
colnames(anno) <- c("Stimulated", "Batch", "Sample", "Material")
rownames(anno) <- colnames(count_mat)

count_mat <- count_mat[order(rownames(count_mat)),]

require("SingleCellExperiment")
sce <- SingleCellExperiment(assays=list(counts=count_mat), colData=anno)
saveRDS(sce, "Tusi_rawCounts.rds")

count_mat <- vector()
anno <- vector();
for (f in files_protein) {
        print(f)
        require("Matrix")
        mat <- read.table(f, sep="\t", header=T)
        mat <- as.matrix(mat)
        active <- "No"
        if (grepl("with_aCD27", f)) {
                active <- "Yes"
        }
        sample <- sub("_matrix.txt.gz","",f)
        sample <- sub("_with_aCD27","",sample)
        sample <- sub("_no_aCD27","",sample)
        sample <- sub("GSM\\d*_[[:alpha:]]*_","",sample)
        biosample <- sub("^\\d*_","",sample)
        a <- cbind(rep(active, ncol(mat)), rep(sample, ncol(mat)), rep(biosample, ncol(mat)), rep("protein", ncol(mat)))
	
	genes <- rownames(mat);
	genes <- unlist(lapply(strsplit(genes, "_"), function(x) {if (length(x) > 2) {return(paste(x, collapse="_"))} else {return(x[1])}}))
	dups <- genes[duplicated(genes)]
	for (d in dups) {
		d_rows <- which(genes == d)
		comb <- colSums(mat[d_rows,])
		mat <- mat[-d_rows,]
		mat <- rbind(mat, comb)
		genes <- genes[-d_rows]
		genes <- c(genes, d);
	}
	rownames(mat) <- genes
	all_g <- sort(union(rownames(count_mat), genes))
	if(!is.null(dim(count_mat))) {
		count_mat <- count_mat[match(all_g, rownames(count_mat)),]
		rownames(count_mat) <- all_g
	}
	mat <- mat[match(all_g, rownames(mat)),]
	rownames(mat) <- all_g

        count_mat <- cbind(count_mat, mat)
        anno <- rbind(anno, a)
}
colnames(anno) <- c("Stimulated", "Batch", "Sample", "Material")
rownames(anno) <- colnames(count_mat)

count_mat <- count_mat[order(rownames(count_mat)),]

controls <- rep(FALSE, nrow(count_mat));
controls[grep("_", rownames(count_mat))] <- TRUE
controls[rownames(count_mat)=="Blank"] <- TRUE

sce_p <- SingleCellExperiment(assays=list(counts=count_mat), colData=anno, rowData=data.frame(feature_symbol=rownames(count_mat), is.control=controls))
saveRDS(sce_p, "Tusi_rawProtein.rds")

ng_thresh = 500
mito_thresh = 0.20
mt_genes <- grep("^MT-", rownames(sce));
ribo_genes <- c(grep("^RPS", rownames(sce)), grep("^RP-", rownames(sce)), grep("^RPL", rownames(sce)))

keepc <- Matrix::colSums(sce@assays[["counts"]] > 0) > ng_thresh & Matrix::colSums(sce@assays[["counts"]][mt_genes,])/Matrix::colSums(sce@assays[["counts"]]) < mito_thresh

sce_clean <- sce[,keepc]
sce_clean <- sce_clean[Matrix::rowSums(sce_clean@assays[["counts"]]) > 10,]
sce_p_clean <- sce_p[,match(colnames(sce_clean), colnames(sce_p))]
sce_p_clean <- sce_p_clean[rowData(sce_p_clean)$is.control==FALSE,]

# hvg = mean > 0.1 & sd > 2 (log-norm)
cannon_markers <- rbind(c("CD3", "T"), c("CD4", "T"), c("CD8", "T"), c("CD11b", "Monocyte"), c("CD33", "Monocyte"), c("CD14", "Monocyte"), c("CD155", "Monocyte"), c("CD19", "B"), c("CD20", "B"), c("CD56", "NK"), c("CD158e1", "NK"), c("PVR", "Monocyte"), c("MS4A1", "B"), c("NCAM1", "NK"), c("KIR3DL1", "NK"), c("FCGR3A", "Monocyte"), c("PF4", "Meg"))

markers <- read.table("/nfs/users/nfs_t/ta6/Data/Hematopoiesis_Markers.txt", header=TRUE)
markers <- markers[order(markers[,1]),]
source("~/R-Scripts/Ensembl_Stuff.R")
markers$Hu_Id <- General_Map(markers[,1], in.org="Mmus", in.name="symbol", out.org="Hsap", out.name="symbol")
markers <- cbind(as.character(markers[,3]), as.character(markers[,2]))
markers <- rbind(markers, cannon_markers)

require("scater")
sce_clean <- scater::normalize(sce_clean)
sce <- sce_clean

g_mean <- Matrix::rowMeans(sce_clean@assays[["logcounts"]])
g_var <- Matrix::rowMeans((as.matrix(sce_clean@assays[["logcounts"]])-g_mean)^2)
sce_clean <- sce_clean[(g_mean > 0.1 & sqrt(g_var)/g_mean > 2) | rownames(sce_clean) %in% cannon_markers[,1] | rownames(sce_clean) %in% markers[,1],]

set.seed(9108)
# 50 PCs - irlba
require(irlba)
pcs <- irlba(sce_clean@assays[["logcounts"]], nv=10)

# k = 40 - Rphenograph
require(Rphenograph)
clusters <- Rphenograph::Rphenograph(pcs$v, 40)

#sce_clean <- sce_clean[rownames(sce_clean) %in% markers[,1],]
#clusters <- kmeans(as.matrix(t(sce_clean@assays[["logcounts"]])), centers=10)


row_mean_aggregate <- function (mat, groups) {
    MAT <- as.matrix(mat)
    x <- split(seq(ncol(MAT)), groups)
    result <- sapply(x, function(a) Matrix::rowMeans(MAT[, a]))
    return(result)
}

profiles <- row_mean_aggregate(sce_clean@assays[["logcounts"]], clusters[[2]]$membership)
#profiles <- row_mean_aggregate(sce_clean@assays[["logcounts"]], clusters$cluster)
rownames(profiles) <- rownames(sce_clean)
profiles <- t(apply(profiles, 1, scale))

tmp <- vector();
for (i in as.character(unique(markers[,2]))) {
	if (sum(rownames(profiles) %in% as.character(markers[markers[,2]==i,1])) > 1) {
        score<- colMeans(profiles[rownames(profiles) %in% as.character(markers[markers[,2]==i,1]),])
	tmp <- rbind(tmp, score)
	} else {
	tmp <- rbind(tmp, rep(0, ncol(tmp)))
	}
}
rownames(tmp) <- as.character(unique(markers[,2]))
tmp[tmp < 0.5] <- 0

groups <- as.character( clusters[[2]]$membership )
for (i in 1:ncol(tmp)) {
	if (sum(tmp[,i]) ==0) {next;}
	g <- which(tmp[,i] == max(tmp[,i]))
	groups[groups == i] <- rownames(tmp)[g]
}

#groups[groups %in% as.character( clusters[[2]]$membership )] <- paste("Undiff", groups[groups %in% as.character( clusters[[2]]$membership )], sep="")

sce$type <- groups;
saveRDS(sce, "Peterson_cleanRNA.rds")
