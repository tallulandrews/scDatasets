# Foreach dataset create an object containing 
#	the QCed expression matrix with unique cell IDs and gene symbols
#	and a cell annotation dataframe with 
#		batch, source, dataset, 
#		known cell-type, and genotype information.

GSE59114_1 = read.delim("GSE59114_C57BL6_GEO_all.csv", header=T) # This file required editing in vi prior to use

expr_mat1 = GSE59114_1[,2:length(GSE59114_1[1,])]
gene_names = as.character(GSE59114_1[,1]);
dups = duplicated(gene_names)==1
gene_names[dups] = paste(gene_names[dups], "1", sep="__")
dups = duplicated(gene_names)==1
gene_names[dups] = paste(gene_names[dups], "1", sep="__")
rownames(expr_mat1) <- gene_names;

keep = colSums(expr_mat1 > 0)>2500;
expr_mat1 = expr_mat1[,keep]
ncells = length(expr_mat1[1,])

labels = as.character(colnames(expr_mat1));

age = rep("old", times=ncells)
age[grep("young",labels)] = "young"
type = rep("Unknown", times=ncells);
type[grep("MPP", labels)] = "Lin-Sca1+Kit+CD150+CD48+"
type[grep("ST_HSC", labels)] = "Lin-Sca1+Kit+CD150-CD48-"
type[grep("LT_HSC", labels)] = "Lin-Sca1+Kit+CD150+CD48-"

colnames(expr_mat1) = paste("GSE59114_C57BL6",colnames(expr_mat1), sep="_")

batch = rep("Unknown", times=ncells);

Annotation1 = data.frame(Cell = colnames(expr_mat1), Batch=batch, Source = rep("MmusBoneMarrow", times=ncells), Dataset = rep("GSE59114_C57BL6", times=ncells), Type = paste(age,type, sep="_"), Genotype=rep("WildType", times=ncells));

Obj = list(data = expr_mat1, anno = Annotation1)
saveRDS(Obj, file="GSE59114_C57BL6.rds")



GSE59114_2 = read.delim("GSE59114_DBA_GEO_all.csv", header=T)
expr_mat2 = GSE59114_2[,2:length(GSE59114_2[1,])]

gene_names = as.character(GSE59114_2[,1]);
dups = duplicated(gene_names)==1
gene_names[dups] = paste(gene_names[dups], "1", sep="__")
rownames(expr_mat2) <- gene_names;

detected = colSums(expr_mat2 > 0)
keep = detected>2500 & detected < 10000;
expr_mat2 = expr_mat2[,keep]
ncells = length(expr_mat2[1,])

labels = as.character(colnames(expr_mat2));

age = rep("Unknown", times=ncells)
age[grep("young",labels)] = "young"
age[grep("old",labels)] = "old"
type = rep("Unknown", times=ncells);
type[grep("MPP", labels)] = "Lin-Sca1+Kit+CD150+CD48+"
type[grep("STHSC", labels)] = "Lin-Sca1+Kit+CD150-CD48-"
type[grep("LTHSC", labels)] = "Lin-Sca1+Kit+CD150+CD48-"

colnames(expr_mat2) = paste("GSE59114_DBA",colnames(expr_mat2), sep="_")

batch = rep("Unknown", times=ncells);

Annotation2 = data.frame(Cell = colnames(expr_mat2), Batch=batch, Source = rep("MmusBoneMarrow", times=ncells), Dataset = rep("GSE59114_DBA", times=ncells), Type = paste(age,type, sep="_"), Genotype=rep("WildType", times=ncells));

Obj = list(data = expr_mat2, anno = Annotation2)
saveRDS(Obj, file="GSE59114_DBA.rds")
