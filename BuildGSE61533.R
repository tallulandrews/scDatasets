# Foreach dataset create an object containing 
#	the QCed expression matrix with unique cell IDs and gene symbols
#	and a cell annotation dataframe with 
#		batch, source, dataset, 
#		known cell-type, and genotype information.

GSE61533 = read.delim("GSE61533_HTSEQ_count_results.csv", header=T)

expr_mat = GSE61533[,2:length(GSE61533[1,])]
rownames(expr_mat) <- GSE61533[,1];

keep = colSums(expr_mat > 0)>4000;
expr_mat = expr_mat[,keep]
ncells = length(expr_mat[1,])
batch = rep("GSE61533",times=ncells)
type = rep("Lin-Kit+Sca1+CD34-Flt3-CD48-CD150+",times=ncells)

Annotation = data.frame(Cell = colnames(expr_mat), Batch=batch, Source = rep("MmusBoneMarrow", times=ncells), Dataset = rep("GSE61533", times=ncells), Type = type, Genotype=rep("WildType", times=ncells));

Obj = list(data = expr_mat, anno = Annotation)
saveRDS(Obj, file="GSE61533.rds")
