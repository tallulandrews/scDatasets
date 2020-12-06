ERCC = read.delim("Zeisel2015ERCC.txt", header=F)
vals = ERCC[-(1:12),-2]
rownames(vals) = vals[,1]
vals = vals[,-1]
ids = as.character(unlist(ERCC[8,-c(1,2)]))
colnames(vals) = ids
ERCC = vals

mtRNA = read.delim("Zeisel2015mtRNA.txt", header=F)
vals = mtRNA[-(1:12),-2]
rownames(vals) = vals[,1]
vals = vals[,-1]
ids = as.character(unlist(mtRNA[8,-c(1,2)]))
colnames(vals) = ids
mtRNA = vals

mRNA = read.delim("Zeisel2015mRNA.txt", header=F)
vals = mRNA[-(1:12),-2]
rownames(vals) = vals[,1]
vals = vals[,-1]
ids = as.character(unlist(mRNA[8,-c(1,2)]))
labels1 = as.character(unlist(mRNA[9,-c(1,2)]))
labels2 = as.character(unlist(mRNA[10,-c(1,2)]))
colnames(vals) = ids
mRNA = vals


mtRNA = mtRNA[,match(colnames(mRNA),colnames(mtRNA))]
data = rbind(mRNA,mtRNA,ERCC)
Zeisel = list(data=data, labels1 = labels1, labels2 = labels2)
