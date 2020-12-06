# Foreach dataset create an object containing 
#	the QCed expression matrix with unique cell IDs and gene symbols
#	and a cell annotation dataframe with 
#		batch, source, dataset, 
#		known cell-type, and genotype information.


####### IN PROGRESS ############

GSE54574 = read.table("GSE54574_E10_E12_HSPC_SingleCell_FPKM.txt", header=T)

expr_mat = GSE54574
keep = colSums(expr_mat > 0) > 2000;
expr_mat = expr_mat[,keep]

labels = colnames(expr_mat)
labels = strsplit(labels, "_")
labels = unlist(lapply(labels, function(x){x[1]}));

ncells = length(expr_mat[1,])
colnames(expr_mat) <- paste("P_",1:ncells, sep="");

batch = paste("GSE54574",labels, sep="_")

type = labels; 
type[type=="HSPC"] = "CD45+Kit+CD34+"; 
type[type=="E10"] = "Prom1+Sca1+CD34+CD45-";
type[type=="E12"] = "Prom1+Sca1+CD34+CD45-";


Annotation = data.frame(Cell = colnames(expr_mat), Batch=batch, Source = rep("MmusPlacenta", times=ncells), Dataset = rep("GSE54574", times=ncells), Type = type, Genotype=rep("WildType", times=ncells));

Obj = list(data = expr_mat, anno = Annotation)
saveRDS(Obj, file="GSE54574.rds")
