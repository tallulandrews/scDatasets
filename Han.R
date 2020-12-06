# Han X, Chen H, Huang D, Chen H et al. Mapping human pluripotent stem cell differentiation pathways using high throughput single-cell RNA-sequencing. Genome Biol 2018 Apr 5;19(1):47.
# PMID : 29622030
# GEO: GSE107552

files <- c("GSM3015987_H20180123-naive.csv",
	"GSM2871123_Primed.csv",
	"GSM2871128_Naive_day10.csv", 
	"GSM2871129_Naive_day20-1.csv", 
	"GSM2871130_Naive_day20-2.csv", 
	"GSM3015985_H20180119.csv", 
	"GSM3015986_H20180123-primed.csv") 
#	"GSM2871124_H20160717_cell_Info.csv", 
#	"GSM2871125_H20160720_cell_info.csv", 
#	"GSM2871126_H20161219_cell_Info.csv", 
#	"GSM2871127_H20161226_cell_Info.csv", 
#	"GSM3015985_H20180119_cell_Info.csv", 
#	"GSM3015986_H20180123_cell_info-primed_H1.csv", 
#	"GSM3015987_H20180123_cell_info-naive_H1.csv"

all <- list()
type <- c();
require("Matrix")
for (f in files) {
	dat <- read.delim(f, sep=",")
	rownames(dat) <- dat[,1]
	dat <- dat[,-1]
	dat <- as.matrix(dat);
	dat <- Matrix(dat);
	type <- c(type, rep(f, times=ncol(dat)))
	all[[f]] <- dat
}

mat <- do.call("cbind", all)
require("SingleCellExperiment")
cell_names <- colnames(mat)
colnames(mat) <- paste("cell", 1:ncol(mat), sep="_");

gene_names <- rownames(mat)
gene_names <- sub("([0-9]+)-Mar", "MARCH\\1", gene_names, perl=TRUE)
gene_names <- sub("([0-9]+)-Dec", "DEC\\1", gene_names, perl=TRUE)
gene_names <- sub("([0-9]+)-Sep", "SEPT\\1", gene_names, perl=TRUE)

sce <- SingleCellExperiment(assay=list(counts=mat), colData=data.frame(file=type, orig_name=cell_names), rowData=data.frame(feature_symbol=rownames(mat)));
saveRDS(sce, "Han.rds")

