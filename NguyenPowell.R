# Nguyen et al. (2018) Single-cell RNA-seq of human induced pluripotent stem cells reveals cellular heterogeneity and cell state transitions between subpopulations. Genome Research. 28:1053-1066
# PMID:29752298
# ArrayExpress: E-MTAB-6687

files = c("hIPSC_scRNA_Sample1.tsv", 
	"hIPSC_scRNA_Sample2.tsv", 
	"hIPSC_scRNA_Sample3.tsv",
	"hIPSC_scRNA_Sample4.tsv",
	"hIPSC_scRNA_Sample5.tsv")


require("Matrix")
data_list <- list()
for( i in 1:length(files) ) {
	f <- files[i]
	dat <- read.table(f, header=T)
	rownames(dat) <- dat[,1]
	dat <- dat[,-1]
	dat <- Matrix(as.matrix(dat))
	cellIDs <- paste("Sample", i, "cell", 1:ncol(dat), sep="_")
	geneids <- rownames(dat)
	gene_names <- matrix(unlist(strsplit(geneids, "_E")), ncol=2, byrow=TRUE)
	ensg <- paste("E", gene_names[,2], sep="")
	rownames(dat) <- ensg;
	colnames(dat) <- cellIDs;
	data_list[[i]] <- dat
}

all <- do.call("cbind", data_list)
require("SingleCellExperiment")
source("~/R-Scripts/Ensembl_Stuff.R")

gene_ids <- data.frame(ensg=rownames(data_list[[1]]), feature_symbol=General_Map(rownames(data_list[[1]]), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol"))

sce <- SingleCellExperiment(assays=list(all), rowData=gene_ids);
saveRDS(sce, file="NguyenPowell.rds")
