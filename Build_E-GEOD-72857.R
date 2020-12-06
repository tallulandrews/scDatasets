# Foreach dataset create an object containing 
#	the QCed expression matrix with unique cell IDs and gene symbols
#	and a cell annotation dataframe with 
#		batch, source, dataset, 
#		known cell-type, and genotype information.

# This is a monster dataset and needs a high-memory job!

source("/nfs/users/nfs_t/ta6/R-Scripts/Name_Map.R")

totcells =  0

DATA = read.table("tmp/GSE72857_umitab.txt", header=T)
Anno = read.delim("tmp/GSE72857_experimental_design.txt", header=T, sep="\t")

expr_mat = as.matrix(DATA);
keep = colSums(expr_mat> 0) > 500;
expr_mat = expr_mat[,keep];

gene_names = rownames(DATA); # What to do with those with many names?
ok = which(!grepl(";", rownames(DATA)))
problem = grep(";", rownames(DATA))

destruct_names <- function(x) {
	y = rownames(DATA)[x]
	possible = unlist(strsplit(y,";"));
	valid = possible[possible %in% map[,2]]
	if (length(valid) == 1) {return(c(1,valid))}
	if (length(valid) < 1) {return(c(1,y))};
	if (length(valid) > 1) {return(c(length(valid),valid))}
}

out = lapply(problem, destruct_names)
reptimes = as.numeric(unlist(lapply(out, function(x) {x[1]})))
rowids = rep(problem, times=reptimes)
my_gene_names = unlist(lapply(out, function(x) {x[-1]}))
expr_mat = expr_mat[c(ok,rowids),]
rownames(expr_mat) <- c(rownames(DATA)[ok],my_gene_names);
## Dealt with: if more than one recognized id in the list duplicate the row for each id

ncells = length(expr_mat[1,])
Anno = Anno[match(colnames(DATA), Anno[,1]),]
Anno = Anno[keep,]
type = rep("Unknown", times=ncells);
geno = rep("WildType", times = ncells); 
type[Anno[,7] =="CMP CD41"] = "CMP";
type[Anno[,7] == "CMP Flt3+ Csf1r+"] = "CMP Flt3+Csf1r+";
type[Anno[,7] == "CMP Irf8-GFP+ MHCII+"] = "CMP Irf8-GFP+MHCII+";
type[Anno[,7] == "Unsorted myeloid"] = "Myeloid";
geno[Anno[,7] =="Cebpa KO"] = "Cebpa KO";
geno[Anno[,7] =="Cebpe KO"] = "Cebpe KO";

batch = Anno[,3]; batch = paste("GEOD72857",batch,sep="_");

Annotation = data.frame(Cell = colnames(expr_mat), Batch=batch, Source = rep("MmusBoneMarrow", times=ncells), Dataset = rep("GEOD76983", times=ncells), Type = type, Genotype=geno);

Obj = list(data = expr_mat, anno = Annotation)
saveRDS(Obj, file="E-GEOD-72857.rds")
