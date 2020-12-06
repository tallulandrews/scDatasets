# Foreach dataset create an object containing 
#	the QCed expression matrix with unique cell IDs and gene symbols
#	and a cell annotation dataframe with 
#		batch, source, dataset, 
#		known cell-type, and genotype information.

GSE60781 = read.delim("GSE60781_RPKM.txt", header=T)

expr_mat = GSE60781[,2:length(GSE60781[1,])]
rownames(expr_mat) <- GSE60781[,1];
type = c(rep("CDP", times=96),rep("PreDC", times=96),rep("MDP",times=59))
batch = paste("GSE60781",type, sep="_")

keep = colSums(expr_mat > 0)>2000;
expr_mat = expr_mat[,keep]
type = type[keep]
batch = batch[keep]
ncells = length(expr_mat[1,])


Annotation = data.frame(Cell = colnames(expr_mat), Batch=batch, Source = rep("MmusBoneMarrow", times=ncells), Dataset = rep("GSE60781", times=ncells), Type = type, Genotype=rep("WildType", times=ncells));

Obj = list(data = expr_mat, anno = Annotation)
saveRDS(Obj, file="GSE60781.rds")
