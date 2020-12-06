# Foreach dataset create an object containing 
#	the QCed expression matrix with unique cell IDs and gene symbols
#	and a cell annotation dataframe with 
#		batch, source, dataset, 
#		known cell-type, and genotype information.

totcells =  0

DATA = read.table("tmp/GSE76983_expdata_BMJhscC.csv", header=T)

expr_mat = as.matrix(DATA[,2:length(DATA[1,])]);
keep = colSums(expr_mat> 0) > 100;
expr_mat = expr_mat[,keep];

gene_names = DATA[,1];
gene_names = matrix(unlist(strsplit(as.character(gene_names),"__")),ncol=2, byrow=T )
dups = duplicated(gene_names[,1])
gene_names[dups,1] = DATA[dups,1]
rownames(expr_mat) <- gene_names[,1];

labels = colnames(expr_mat);
ncells = length(labels)
type = rep("Unknown", times=ncells); type[grepl("HSC",labels)] <- "Lin-Sca1+Kit+CD150+CD48-"; #HSC
thing = strsplit(labels,"_")
batch = unlist(lapply(thing, function(x){x[1]})); batch = paste("GEOD76983",batch,sep="_");

Annotation = data.frame(Cell = colnames(expr_mat), Batch=batch, Source = rep("MmusBoneMarrow", times=ncells), Dataset = rep("GEOD76983", times=ncells), Type = type, Genotype=rep("WildType", times=ncells));

Obj = list(data = expr_mat, anno = Annotation)
saveRDS(Obj, file="E-GEOD-76983.rds")
