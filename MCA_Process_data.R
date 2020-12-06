# R3.4.2
args <- commandArgs(trailingOnly=TRUE) ### For Farm

require("Matrix")
require("SingleCellExperiment")

gene_rule="intersect" # rule for homogenizing genes

ann <- read.delim("MCA_All-batch-removed-assignments.csv", header=T)
tissues <- lapply(strsplit(as.character(ann[,1]), "_"), function(x){x[1]})
tissues <- unique(unlist(tissues))
tissues <- tissues[tissues != "BoneMarrowcKit"]
tissues[tissues=="Male(fetal)Gonad"] <- "FetalMaleGonad"
tissues[tissues=="NeonatalBrain"] <- "NeontalBrain"

#tissues <- levels(ann$Tissue)
#tissues <- gsub("-","", tissues)
#tissues <- gsub("_","", tissues)

#all_files <- Sys.glob(paste("MCA_DGE/*.txt.gz", sep=""))
#for(t in tissues) {
#	exclude <- grep(t, all_files)
#	if (length(exclude) > 0) {
#		all_files <- all_files[-exclude];
#	}
#}
#extras <- sub("_dge.txt.gz", "", all_files)
#extras <- sub("MCA_DGE/", "", extras)
#tissues <- c(tissues, extras)

# "NeontalBrain"

tissues <- tissues[as.numeric(args[1])] ### For Farm

for (t in tissues) {
	files1 <- Sys.glob(paste("MCA_DGE/",t,"*.txt.gz", sep="")) 
	files2 <- Sys.glob(paste("MCA_rmbatch_dge_all/",t,"*.txt.gz",sep=""))

	if (length(files1) < 1) {next;}
	sparsemats <- list()
	annmats <-list()
	for (f in files1) {
		dat <- read.delim(f, sep=" ", header=T)
		M <- Matrix::Matrix(as.matrix(dat), sparse=T)
		local_ann <- ann[match(colnames(dat), ann[,1]),]
		local_ann$Cell.name <- colnames(dat)
		local_ann$ClusterID <- as.character(local_ann$ClusterID)
		local_ann$ClusterID[is.na(local_ann$ClusterID)] <- "unassigned"
		local_ann$ClusterID<- factor(local_ann$ClusterID)
		local_ann$Tissue[is.na(local_ann$Tissue)] <- local_ann$Tissue[which(!is.na(local_ann$Tissue))[1]]
		local_ann$Batch[is.na(local_ann$Batch)] <- NA
		local_ann <- local_ann[,-ncol(local_ann)]
		rownames(local_ann) <- local_ann[,1]
		local_ann <- local_ann[,-1]
		local_ann$file <- rep(f, times=nrow(local_ann))

		sparsemats[[f]] <- M
		annmats[[f]] <- local_ann
	}
	names(annmats) <- rep(NULL, times=length(annmats))
	names(sparsemats) <- rep(NULL, times=length(sparsemats))

	batchmats <- list()
	for (f in files2) {
		dat <- read.delim(f, sep=" ", header=T)
                M <- Matrix::Matrix(as.matrix(dat), sparse=T)
		batchmats[[f]] <- M
	}

	# Homogenize genes across files for same tissue
	genes <- c();
	for (i in c(sparsemats, batchmats)) {
		if (length(genes) == 0) {
			genes <- rownames(i)
		} else {
			if (gene_rule == "intersect") {
				genes <- intersect(genes, rownames(i))
			} else if (gene_rule == "union") {
				genes <- union(genes, rownames(i))
			} else {
				warn("I don't know how to homogenize genes - assuming all have the same")
				break;
			}
		}
	}
	genes <- sort(genes);

	for(i in 1:length(sparsemats)) {
		matches <- match(genes, rownames(sparsemats[[i]]))
		sparsemats[[i]] <- rbind( sparsemats[[i]], rep(0, times=ncol(sparsemats[[i]])) )
		matches[is.na(matches)] <- nrow(sparsemats[[i]]);

		sparsemats[[i]] <- sparsemats[[i]][ matches, ]
		rownames(sparsemats[[i]]) <- genes
	}
	for(i in 1:length(batchmats)) {
		matches <- match(genes, rownames(batchmats[[i]]))
		batchmats[[i]] <- rbind( batchmats[[i]], rep(0, times=ncol(batchmats[[i]])) )
		matches[is.na(matches)] <- nrow(batchmats[[i]]);

		batchmats[[i]] <- batchmats[[i]][ matches, ]
		rownames(batchmats[[i]]) <- genes
	}

	all_molecules <- do.call("cbind", sparsemats)
	all_batchcor <- do.call("cbind", batchmats)
	all_anns <- do.call("rbind", annmats)
	if (!identical(rownames(all_anns), colnames(all_molecules))) {print("Something is wrong")}
	if (!identical(colnames(all_batchcor), colnames(all_molecules))) {
		matches <- match(colnames(all_molecules), colnames(all_batchcor))
		all_batchcor <- cbind(all_batchcor, rep(0, times=nrow(all_batchcor)))
		matches[is.na(matches)] <- ncol(all_batchcor)

		all_batchcor <- all_batchcor[,matches]
		colnames(all_batchcor) <- colnames(all_molecules)
	}


	sceset <- SingleCellExperiment(assays = list(counts = all_molecules), colData=all_anns)
	assay(sceset, "batchcor") <- all_batchcor

	print(paste(t, "Done!"))
	saveRDS(sceset, file=paste(t,"SingCellExp_sparseM.rds", sep="_"))
}
