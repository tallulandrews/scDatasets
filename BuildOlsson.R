# Foreach dataset create an object containing 
#	the QCed expression matrix with unique cell IDs and gene symbols
#	and a cell annotation dataframe with 
#		batch, source, dataset, 
#		known cell-type, and genotype information.

totcells =  0

Olsson1 = read.table("Olsson_BM_lineageNeg_CD117pos_CD34pos.txt", header=T)
expr_mat = Olsson1[,2:length(Olsson1[1,])];
rownames(expr_mat) = Olsson1[,1];
expr_mat = expr_mat[,colSums(expr_mat > 0)>2000]
colnames(expr_mat) = paste("O_BM_", 1:length(expr_mat[1,]), sep="")
Olsson1 <- expr_mat;
ncells = length(Olsson1[1,])
totcells =  totcells+ncells

Annotation1 = data.frame(Cell = colnames(Olsson1), Batch=rep("Olsson1", times=ncells), Source = rep("MmusBoneMarrow", times=ncells), Dataset = rep("Olsson", times=ncells), Type = rep("Lin-CD117+CD34+", times=ncells), Genotype=rep("WildType", times=ncells));

Olsson2 = read.table("Olsson_BM_lineageNeg_CD117pos_Sca1pos.txt", header=T)
expr_mat = Olsson2[,2:length(Olsson2[1,])];
rownames(expr_mat) = Olsson2[,1];
expr_mat = expr_mat[,colSums(expr_mat > 0)>2000]
colnames(expr_mat) = paste("O_BM_", totcells+1:length(expr_mat[1,]), sep="");
Olsson2 <- expr_mat;
ncells = length(Olsson2[1,]);
totcells =  totcells+ncells

Annotation2 = data.frame(Cell = colnames(Olsson2), Batch=rep("Olsson2", times=ncells), Source = rep("MmusBoneMarrow", times=ncells), Dataset = rep("Olsson", times=ncells), Type = rep("Lin-CD117+Sca1+", times=ncells), Genotype=rep("WildType", times=ncells));


Olsson3 = read.table("Olsson_CMP.txt", header=T)
expr_mat = Olsson3[,2:length(Olsson3[1,])];
rownames(expr_mat) = Olsson3[,1];
expr_mat = expr_mat[,colSums(expr_mat > 0)>2000]
colnames(expr_mat) = paste("O_CMP_", totcells+1:length(expr_mat[1,]), sep="");
Olsson3 <- expr_mat;
ncells = length(Olsson3[1,]);
totcells =  totcells+ncells

Annotation3 = data.frame(Cell = colnames(Olsson3), Batch=rep("Olsson3", times=ncells), Source = rep("MmusBoneMarrow", times=ncells), Dataset = rep("Olsson", times=ncells), Type = rep("CMP", times=ncells), Genotype=rep("WildType", times=ncells));


Olsson4 = read.table("Olsson_GMP.txt", header=T)
expr_mat = Olsson4[,2:length(Olsson4[1,])];
rownames(expr_mat) = Olsson4[,1];
expr_mat = expr_mat[,colSums(expr_mat > 0)>2000]
colnames(expr_mat) = paste("O_GMP_", totcells+1:length(expr_mat[1,]), sep="");
Olsson4 <- expr_mat;
ncells = length(Olsson4[1,]);
totcells =  totcells+ncells

Annotation4 = data.frame(Cell = colnames(Olsson4), Batch=rep("Olsson4", times=ncells), Source = rep("MmusBoneMarrow", times=ncells), Dataset = rep("Olsson", times=ncells), Type = rep("GMP", times=ncells), Genotype=rep("WildType", times=ncells));


Olsson5 = read.table("Olsson_Gfi1-GFP_GMP.txt", header=T)
expr_mat = Olsson5[,2:length(Olsson5[1,])];
rownames(expr_mat) = Olsson5[,1];
expr_mat = expr_mat[,colSums(expr_mat > 0)>2000]
colnames(expr_mat) = paste("O_GMP_", totcells+1:length(expr_mat[1,]), sep="");
Olsson5 <- expr_mat;
ncells = length(Olsson5[1,]);
totcells =  totcells+ncells

Annotation5 = data.frame(Cell = colnames(Olsson5), Batch=rep("Olsson5", times=ncells), Source = rep("MmusBoneMarrow", times=ncells), Dataset = rep("Olsson", times=ncells), Type = rep("GMP", times=ncells), Genotype=rep("Gfi1-GFP", times=ncells));


Olsson6 = read.table("Olsson_Irf8-GFP_GMP.txt", header=T)
expr_mat = Olsson6[,2:length(Olsson6[1,])];
rownames(expr_mat) = Olsson6[,1];
expr_mat = expr_mat[,colSums(expr_mat > 0)>2000]
colnames(expr_mat) = paste("O_GMP_", totcells+1:length(expr_mat[1,]), sep="");
Olsson6 <- expr_mat;
ncells = length(Olsson6[1,]);
totcells =  totcells+ncells

Annotation6 = data.frame(Cell = colnames(Olsson6), Batch=rep("Olsson6", times=ncells), Source = rep("MmusBoneMarrow", times=ncells), Dataset = rep("Olsson", times=ncells), Type = rep("GMP", times=ncells), Genotype=rep("Irf8-GFP", times=ncells));


Olsson7 = read.table("Olsson_Irf8.Null.txt", header=T)
expr_mat = Olsson7[,2:length(Olsson7[1,])];
rownames(expr_mat) = Olsson7[,1];
expr_mat = expr_mat[,colSums(expr_mat > 0)>2000]
colnames(expr_mat) = paste("O_Irf8_", totcells+1:length(expr_mat[1,]), sep="");
Olsson7 <- expr_mat;
ncells = length(Olsson7[1,]);
totcells =  totcells+ncells

Annotation7 = data.frame(Cell = colnames(Olsson7), Batch=rep("Olsson7", times=ncells), Source = rep("MmusBoneMarrow", times=ncells), Dataset = rep("Olsson", times=ncells), Type = rep("Unknown", times=ncells), Genotype=rep("Irf8-KO", times=ncells));


Olsson8 = read.table("Olsson_Gfi1.Null.txt", header=T)
expr_mat = Olsson8[,2:length(Olsson8[1,])];
rownames(expr_mat) = Olsson8[,1];
expr_mat = expr_mat[,colSums(expr_mat > 0)>2000]
colnames(expr_mat) = paste("O_Gfi1_", totcells+1:length(expr_mat[1,]), sep="");
Olsson8 <- expr_mat;
ncells = length(Olsson8[1,]);
totcells =  totcells+ncells

Annotation8 = data.frame(Cell = colnames(Olsson8), Batch=rep("Olsson8", times=ncells), Source = rep("MmusBoneMarrow", times=ncells), Dataset = rep("Olsson", times=ncells), Type = rep("Unknown", times=ncells), Genotype=rep("Gfi1-KO", times=ncells));


Olsson9 = read.table("Olsson_Gfi1.Null.Irf8.Null.txt", header=T)
expr_mat = Olsson9[,2:length(Olsson9[1,])];
rownames(expr_mat) = Olsson9[,1];
expr_mat = expr_mat[,colSums(expr_mat > 0)>2000]
colnames(expr_mat) = paste("O_Irf8_Gfi1_", totcells+1:length(expr_mat[1,]), sep="");
Olsson9 <- expr_mat;
ncells = length(Olsson9[1,]);
totcells =  totcells+ncells

Annotation9 = data.frame(Cell = colnames(Olsson9), Batch=rep("Olsson9", times=ncells), Source = rep("MmusBoneMarrow", times=ncells), Dataset = rep("Olsson", times=ncells), Type = rep("Unknown", times=ncells), Genotype=rep("Gfi1-KO,Irf8-KO", times=ncells));

Annotation = rbind(Annotation1,Annotation2,Annotation3,Annotation4,Annotation5, Annotation6, Annotation7, Annotation8, Annotation9);

Data = cbind(Olsson1, Olsson2, Olsson3, Olsson4, Olsson5, Olsson6, Olsson7, Olsson8, Olsson9)

Obj = list(data = Data, anno = Annotation)
saveRDS(Obj, file="Olsson.rds")
