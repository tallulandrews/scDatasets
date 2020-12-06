# Foreach dataset create an object containing 
#	the QCed expression matrix with unique cell IDs and gene symbols
#	and a cell annotation dataframe with 
#		batch, source, dataset, 
#		known cell-type, and genotype information.

totcells =  0

Drissen = read.table("Drissen_et.al_RefSeq_FPKM.txt", header=T)


expr_mat = Drissen[,2:length(Drissen[1,])]
rownames(expr_mat) <- Drissen[,1];
colnames(expr_mat) <- paste("D_",colnames(expr_mat), sep="");

ncells = length(expr_mat[1,])

Annotation = data.frame(Cell = colnames(expr_mat), Batch=rep("Drissen1",times=ncells), Source = rep("MmusBoneMarrow", times=ncells), Dataset = rep("Drissen", times=ncells), Type = rep("pre-GM", times=ncells), Genotype=rep("WildType", times=ncells));

Obj = list(data = expr_mat, anno = Annotation)
saveRDS(Obj, file="Drissen.rds")
