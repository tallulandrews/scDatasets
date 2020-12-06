# Foreach dataset create an object containing 
#	the QCed expression matrix with unique cell IDs and gene symbols
#	and a cell annotation dataframe with 
#		batch, source, dataset, 
#		known cell-type, and genotype information.

GSE67120 = read.delim("GSE67120_Fan_HSPC_fpkms.out", header=T)

expr_mat = GSE67120[,2:length(GSE67120[1,])]
rownames(expr_mat) <- GSE67120[,1];
colnames(expr_mat) <- paste("GSE67120",colnames(expr_mat), sep="_")

keep = colSums(expr_mat > 0)>4000;
expr_mat = expr_mat[,keep]
ncells = length(expr_mat[1,])
batch = rep("Unknown",times=ncells)

labels = colnames(expr_mat)
my_source = rep("Mmus aorta-gonad-mesonephros", times=ncells)
my_source[grep("FL", labels)] = "Mmus fetal liver"
my_source[grep("Adult",labels)] = "MmusBoneMarrow";

type = rep("Unknown", times =  ncells)
type[grep("EC", labels)] = "CD31+VE-cadherin+CD41-CD43-CD45-Ter119-";
type[grep("_T1_", labels)] = "CD31+CD45-CD41loKit+CD201hi";
type[grep("_T2_", labels)] = "CD31+CD45+Kit+CD201hi"; # I think this is right.
type[grep("_T1CD201neg", labels)] = "CD31+CD45-CD41loKit+CD201-";
type[grep("E14", labels)] = "CD45+CD150+CD48-CD201+";
type[grep("E12", labels)] = "Lin-Sca1+Mac-1loCD201+";
type[grep("Adult", labels)] = "AdultHSC";

# aorta-gonad-mesonephros region (AGM) = endothelial, pre-HSCs
# liver = HSCs

# pre-HSC = CD31+CD45-CD41lowKit+CD201hi

#AGM
# ECs = CD31+VE-cadherin+CD41-CD43-CD45-Ter119-
# T1 = CD31+CD45-CD41loKit+CD201hi
# T2 = CD31+CD45+CD41lo (may also contain CD31+CD45+Kit+CD201hi - or these may be in the ECs?)

#liver
# E12 HSCs = Lin-Sca1+Mac-1loCD201+
# E14 HSCs = CD45+CD150+CD48-CD201+



Annotation = data.frame(Cell = colnames(expr_mat), Batch=batch, Source = rep("MmusBoneMarrow", times=ncells), Dataset = rep("GSE67120", times=ncells), Type = type, Genotype=rep("WildType", times=ncells));

Obj = list(data = expr_mat, anno = Annotation)
saveRDS(Obj, file="GSE67120.rds")
