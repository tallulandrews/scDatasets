# Foreach dataset create an object containing 
#	the QCed expression matrix with unique cell IDs and gene symbols
#	and a cell annotation dataframe with 
#		batch, source, dataset, 
#		known cell-type, and genotype information.

Grover = read.table("Grover.A_et.al_RefSeq.Read.Count.txt", header=T)

expr_mat = Grover
labels = colnames(Grover)
keep = colSums(expr_mat > 0) > 3000;

labels = unlist(strsplit(labels, "_"))
labels = labels[seq(from=2, to=length(labels), by=2)]

expr_mat = expr_mat[,keep]
labels = labels[keep]
ncells = length(expr_mat[1,])
colnames(expr_mat) <- paste("G_",1:ncells, sep="");

batch = rep("Grover1",times=ncells); batch[labels=="old"] <- "Grover2";
type = rep("Lin-Sca1+Kit+CD150+CD48-", times=ncells) #HSC
type[labels == "old"] = "old_Lin-Sca1+Kit+CD150+CD48-";


Annotation = data.frame(Cell = colnames(expr_mat), Batch=batch, Source = rep("MmusBoneMarrow", times=ncells), Dataset = rep("Grover", times=ncells), Type = type, Genotype=rep("WildType", times=ncells));

Obj = list(data = expr_mat, anno = Annotation)
saveRDS(Obj, file="Grover.rds")
