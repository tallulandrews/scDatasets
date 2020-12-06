# Foreach dataset create an object containing 
#	the QCed expression matrix with unique cell IDs and gene symbols
#	and a cell annotation dataframe with 
#		batch, source, dataset, 
#		known cell-type, and genotype information.

#How to convert gene ids?
source("/nfs/users/nfs_t/ta6/R-Scripts/Name_Map.R")

totcells =  0

Nestorowa = read.table("Nestorowa_expression_matrix.txt.gz", header=T)
Anno = read.table("Nestorowa_cell_type_annotation.txt.gz", header=T)


expr_mat = t(Nestorowa[,28:length(Nestorowa[1,])])
rownames(expr_mat) <- ensg2symbol(rownames(expr_mat));
colnames(expr_mat) <- paste("N_",Nestorowa[,27], sep="");
type <- Nestorowa[,26];
batch <- Nestorowa[,25];
batch = as.character(batch);
batch[batch == "batch_1"] <- "Nestorowa1"
batch[batch == "batch_2"] <- "Nestorowa2"

celltype = Anno[,c(13,14,16:22)]
celltype_names = c("LTHSC", "LMPP", "CMP", "MEP","GMP","MPP","MPP","MPP","STHSC")
type2 = apply(celltype, 1, function(x){if(sum(x) == 0) {return("Unknown")} else {return(celltype_names[which(x==1)])}});

type2 = type2[match(Nestorowa[,27], names(type2))]

keep = !is.na(type2) & colSums(expr_mat > 0) > 2000
type2 = type2[keep]
batch = batch[keep]
expr_mat = expr_mat[,keep]
ncells = length(batch)

#Obj = list(data = Data, anno = Annotation)
#saveRDS(Obj, file="Olsson.rds")

Annotation = data.frame(Cell = colnames(expr_mat), Batch=batch, Source = rep("MmusBoneMarrow", times=ncells), Dataset = rep("Nestorowa", times=ncells), Type = type2, Genotype=rep("WildType", times=ncells));
rownames(Annotation) <- colnames(expr_mat)

Obj = list(data = expr_mat, anno = Annotation)
require("SingleCellExperiment")
require("Matrix")
Obj <- SingleCellExperiment(assays=list(counts=Matrix(as.matrix(expr_mat))), colData=Annotation)
saveRDS(Obj, file="Nestorowa.rds")
