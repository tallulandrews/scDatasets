# Estimated purity:
# monocytes: 98%
# B cells: 100%
# CD34+: 45% 
# CD4+ Helper T cells: 99%
# CD4+CD25+ T-reg: 95%
# CD4+CD45RA+/CD25- Naive T cells: 98%
# CD4+/CD45RO+ Memory T cells: 98%
# CD56+ Natural Killer Cells: 92%
# CD8+ Cytotoxic T cells: 98%
# CD8+/CD45RA+ Naive Cytotoxic T cells: 99%

set.seed(134)

arg_list <- c("cd14_monocytes_filtered_gene_bc_matrices.tar.gz",
		"b_cells_filtered_gene_bc_matrices.tar.gz",
#		"cd34_filtered_gene_bc_matrices.tar.gz",
		"cd4_t_helper_filtered_gene_bc_matrices.tar.gz",
		"regulatory_t_filtered_gene_bc_matrices.tar.gz",
		"naive_t_filtered_gene_bc_matrices.tar.gz",
		"memory_t_filtered_gene_bc_matrices.tar.gz",
		"cd56_nk_filtered_gene_bc_matrices.tar.gz",
		"cytotoxic_t_filtered_gene_bc_matrices.tar.gz",
		"naive_cytotoxic_filtered_gene_bc_matrices.tar.gz")
Reference_types <- c("monocyte", "b_cell", "t_helper", "t_reg", "t_naive", "t_mem", "nk", "t_cyto", "t_naive_cyto")

Profiles <- list();
Subset <- list();
for (i in 1:length(arg_list)) {
#args <- commandArgs(trailingOnly=TRUE)
	system(paste("tar -xvzf ",arg_list[i]));

	files <- c("filtered_matrices_mex/hg19/matrix.mtx","filtered_matrices_mex/hg19/barcodes.tsv", "filtered_matrices_mex/hg19/genes.tsv")
	counts <- read.delim(files[1], stringsAsFactors=FALSE, sep=" ", header=FALSE)
	barcode <- read.delim(files[2], stringsAsFactors=FALSE, header=FALSE)
	genes <- read.delim(files[3], stringsAsFactors=FALSE, header=FALSE)
	
	comments <- grep("%", counts[,1])
	counts <- counts[-comments,]
	counts <- counts[!is.na(counts[,1]) & !is.na(counts[,2]) & !is.na(counts[,3]),]
	#counts[,3] <- log2(as.numeric(counts[,3]) +1) # for randomforest


	require("Matrix")
	mat <- sparseMatrix(i=as.numeric(counts[,1]), j=as.numeric(counts[,2]), x = as.numeric(counts[,3]), dimnames = list(genes[,2], barcode[,1]))
	
	Data <- list(counts = mat, type = Reference_types[i])
	parsed_filename <- strsplit(arg_list[i], "filtered")

	saveRDS(Data, file=paste(parsed_filename[[1]][1], "dataset.rds", sep=""))

	Profiles[[i]] <- rowMeans(mat)
	exclude <- unique(sort(c(grep("^RPS",rownames(mat)), grep("^RPL", rownames(mat)), grep("^MT-", rownames(mat)), which(Profiles[[i]] < 0.1), which(rownames(mat)=="AC002321.1"))))

	require("proxy")
	res <- proxy::simil(t(Profiles[[i]][-exclude]), t(as.matrix(mat[-exclude,])), method="cosine")
	keep <- which(res > quantile(res, probs=0.9))
	set.seed(141)
	keep <- runif(n=ncol(mat)) < 0.1; # Keep 10% sample
	Subset[[i]] <- mat[,keep]
}

P_Mat <- cbind(Profiles[[1]],Profiles[[2]], Profiles[[3]], Profiles[[4]], Profiles[[5]], Profiles[[6]], Profiles[[7]], Profiles[[8]], Profiles[[9]])
#colnames(P_Mat) <- c("Mono", "B cells", "T-help", "T-reg", "T-naive", "T-mem", "NK", "T-cyto", "T-naive-cyto")
colnames(P_Mat) <- Reference_types
sf <- colSums(P_Mat)
P_Mat <- t(t(P_Mat)/sf*median(sf))
P_Mat <- log2(P_Mat+1)

P_max <- apply(P_Mat, 1, max)
P_min <- apply(P_Mat, 1, min)

saveRDS(P_Mat, file="reference_logmat.rds")

Features <- P_max/P_min > 2 & P_max > 0.1 
# T-cell Features
T_mat <- P_Mat[,-c(1,2,7)]
Features2 <- apply(T_mat,1,max)/apply(T_mat,1,min) & P_max > 0.1

cosineDist <- function(x){
  as.dist(1 - x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))) 
}
cosineSim <- function(x){
  x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))
}

# Fancy Feature selection (Balanced)
Potential <- P_max > 0.05
Pot_Mat <- P_Mat[Potential,]
score <- P_max[Potential]/P_min[Potential]
ID <- apply(P_Mat, 1, function(x){which(x==max(x))[1]})
ID <- ID[Potential]
reorder <- order(-score)
ID <- ID[reorder]
Features3 <- vector();
for (i in 1:ncol(P_Mat)) {
	Features3 <- c(Features3,names(ID[ID==i])[1:200])
}
Features3 <- rownames(P_Mat) %in% Features3


# Check Assignments
Assigned <- matrix(0, ncol=length(arg_list), nrow=length(arg_list))
for(i in 1:length(arg_list)) {
	 system(paste("tar -xvzf ",arg_list[i]));

        files <- c("filtered_matrices_mex/hg19/matrix.mtx","filtered_matrices_mex/hg19/barcodes.tsv", "filtered_matrices_mex/hg19/genes.tsv")
        counts <- read.delim(files[1], stringsAsFactors=FALSE, sep=" ", header=FALSE)
        barcode <- read.delim(files[2], stringsAsFactors=FALSE, header=FALSE)
        genes <- read.delim(files[3], stringsAsFactors=FALSE, header=FALSE)

        comments <- grep("%", counts[,1])
        counts <- counts[-comments,]
        counts <- counts[!is.na(counts[,1]) & !is.na(counts[,2]) & !is.na(counts[,3]),]

        require("Matrix")
        mat <- sparseMatrix(i=as.numeric(counts[,1]), j=as.numeric(counts[,2]), x = as.numeric(counts[,3]), dimnames = list(genes[,2], barcode[,1]))

	require("proxy")
	res <- proxy::simil(t(P_Mat[Features,]), t(as.matrix(mat[Features,])), method="cosine")
	correct <- res[i,] == apply(res,2,max)
	assign <- apply(res,2,function(x){if(max(x) > 0.7){which(x == max(x))[1]} else{NA}})

	factor_counts <- function(vec) {
		 x <- split(seq(length(vec)), vec)
		 result <- sapply(x, function(a) length(a))
		 return(result)
	}

	Assigned[i,] <- factor_counts(factor(assign, levels=1:9))
}	

# RandomForest

test <- cbind(as.matrix(Subset[[1]]), as.matrix(Subset[[2]]), as.matrix(Subset[[3]]), as.matrix(Subset[[4]]), as.matrix(Subset[[5]]), as.matrix(Subset[[6]]), as.matrix(Subset[[7]]), as.matrix(Subset[[8]]), as.matrix(Subset[[9]]))
test_lab <- rep( c("Mono", "B cells", "T-help", "T-reg", "T-naive", "T-mem", "NK", "T-cyto", "T-naive-cyto"), times=c(ncol(Subset[[1]]), ncol(Subset[[2]]), ncol(Subset[[3]]), ncol(Subset[[4]]), ncol(Subset[[5]]), ncol(Subset[[6]]), ncol(Subset[[7]]), ncol(Subset[[8]]), ncol(Subset[[9]])))
require("randomForest")
res <- randomForest(t(test), factor(test_lab), ntree=50, keep.forest=TRUE) # SLoW
# still 30-60% error rate for T-cell sub-types



Features_names <- rownames(P_Mat)[Features]


### Assign Data using Reference ###

to_assign <- c("fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz", 
		"frozen_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz", 
		"frozen_pbmc_donor_b_filtered_gene_bc_matrices.tar.gz",
		"frozen_pbmc_donor_c_filtered_gene_bc_matrices.tar.gz",
		"pbmc33k_filtered_gene_bc_matrices.tar.gz",
		"pbmc3k_filtered_gene_bc_matrices.tar.gz",
		"pbmc4k_filtered_gene_bc_matrices.tar.gz",
		"pbmc6k_filtered_gene_bc_matrices.tar.gz",
		"pbmc8k_filtered_gene_bc_matrices.tar.gz")
Assigned_new <- matrix(0, ncol=length(arg_list), nrow=length(to_assign))
for (i in 1:length(to_assign)) {
	system(paste("tar -xvzf ", to_assign[i], "> tmp.txt"));
	files <- as.vector(read.table("tmp.txt")[,1])
	mat_file <- files[grep("matrix.mtx", files)]
	bar_file <- files[grep("barcodes", files)]
	gen_file <- files[grep("genes", files)]

	counts <- read.delim(mat_file, stringsAsFactors=FALSE, sep=" ", header=FALSE)
	barcode <- read.delim(bar_file, stringsAsFactors=FALSE, header=FALSE)
	genes <- read.delim(gen_file, stringsAsFactors=FALSE, header=FALSE)

	comments <- grep("%", counts[,1])
        counts <- counts[-comments,]
        counts <- counts[!is.na(counts[,1]) & !is.na(counts[,2]) & !is.na(counts[,3]),]

        require("Matrix")
        mat <- sparseMatrix(i=as.numeric(counts[,1]), j=as.numeric(counts[,2]), x = as.numeric(counts[,3]), dimnames = list(genes[,2], barcode[,1]))
	mat <- mat[rownames(mat) %in% rownames(P_Mat),]
	P_Mat_tmp <- P_Mat[match(rownames(mat), rownames(P_Mat)),]

	# Get Cluster data
	parsed_filename <- strsplit(to_assign[i], "filtered")
	clust_file <- paste(parsed_filename[[1]][1], "analysis_k10.csv.gz", sep="")

	comp_clust <- read.table(clust_file, stringsAsFactors=FALSE, header=TRUE, sep=",")

	# Assign labels by reference
	require("proxy")
	res <- proxy::simil(t(P_Mat_tmp[rownames(P_Mat_tmp) %in% Features_names,]), 
			    t(as.matrix(mat[rownames(mat) %in% Features_names,])), method="cosine")
	correct <- res[i,] == apply(res,2,max)
	assign <- apply(res,2,function(x){if(max(x) > 0.7){which(x == max(x))[1]} else{NA}})

	Data <- list(counts = mat, clusters = comp_clust[,2], assigned = Reference_types[assign])
	saveRDS(Data, file=paste(parsed_filename[[1]][1], "dataset.rds", sep=""))

	Assigned_new[i,] <- factor_counts(factor(assign, levels=1:9))
}


























args <- commandArgs(trailingOnly=TRUE)

args <- c("filtered_gene_bc_matrices/hg19/matrix.mtx","filtered_gene_bc_matrices/hg19/barcodes.tsv","filtered_gene_bc_matrices/hg19/genes.tsv", "pbmc3k_analysis_k10.csv.gz", "out")

counts <- read.delim(args[1], stringsAsFactors=FALSE, sep=" ", header=FALSE)
barcode <- read.delim(args[2], stringsAsFactors=FALSE, header=FALSE)
genes <- read.delim(args[3], stringsAsFactors=FALSE, header=FALSE)
clust <- read.delim(args[4], header=TRUE, sep=",", stringsAsFactors=FALSE)

out <- args[5];

comments <- grep("%", counts[,1])
counts <- counts[-comments,]
counts <- counts[!is.na(counts[,1]) & !is.na(counts[,2]) & !is.na(counts[,3]),]

require("Matrix")
mat <- sparseMatrix(i=as.numeric(counts[,1]), j=as.numeric(counts[,2]), x = as.numeric(counts[,3]), dimnames = list(genes[,1], barcode[,1]))

require("CellTypeProfiles")

# auto_QC 
Ts <- colSums(mat)
Ds <- colSums(mat > 0)
tot <- quantile(Ts, probs=c(0.25, 0.5, 0.75))
fea <- quantile(Ds, probs=c(0.25, 0.5, 0.75))
#outliers1 <- abs(Ts - tot[2]) > 3.5*(mean(tot[2]-tot[1], tot[3]-tot[2]))
#outliers2 <- abs(Ds - fea[2]) > 3.5*(mean(fea[2]-fea[1], fea[3]-fea[2]))
outliers1 <- abs(Ts - tot[2]) > (tot[2]-min(Ts))
outliers2 <- abs(Ds - fea[2]) > (fea[2]-min(Ds))

mat <- mat[, !(outliers1 | outliers2)]
clust <- clust[!(outliers1 | outliers2),]

clust_sizes <- factor_counts(clust[,2])
exclude <- names(clust_sizes)[clust_sizes<=1]
mat <- mat[,!clust[,2] %in% exclude]
clust <- clust[!clust[,2] %in% exclude,]

profiles <- CellTypeProfiles::get_cluster_profile(mat, clust[,2], feature_selection=m3drop.features, norm_method="CPM")

saveRDS(profiles, file=paste(out,"_profile.rds", sep=""))
saveRDS(mat, file=paste(out,"_mat.rds", sep=""))

