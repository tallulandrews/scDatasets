# Kolodziejczyk et al. (2015). Single Cell RNA-Sequencing of Pluripotent States Unlocks Modular Transcriptional Variation. Cell Stem Cell. 17(4) 471-485
# PMID: 26431182

data <- read.table("counttable_es.csv", header=TRUE)

anno = colnames(data);
anno[grep("_2i_", anno)] = "2i"
anno[grep("_a2i_", anno)] = "a2i"
anno[grep("_lif_", anno)] = "lif"

batch = colnames(data);
batch[grep("_2i_1", batch)] = "2i_1"
batch[grep("_2i_2", batch)] = "2i_2"
batch[grep("_2i_3", batch)] = "2i_3"
batch[grep("_a2i_1", batch)] = "a2i_1"
batch[grep("_a2i_2", batch)] = "a2i_2"
batch[grep("_a2i_3", batch)] = "a2i_3"
batch[grep("_lif_1", batch)] = "lif_1"
batch[grep("_lif_2", batch)] = "lif_2"
batch[grep("_lif_3", batch)] = "lif_3"

spikes <- grep("ERCC-", rownames(data), ignore.case=TRUE)




ANN <- data.frame(Species = rep("Mus musculus", times=length(anno)), cell_type1 = anno, batch=batch, Source=rep("ESC", times=length(anno)))
rownames(ANN) <- colnames(data);

require("scater")
pd <- new("AnnotatedDataFrame", data=ANN)
kolo <- newSCESet(countData=as.matrix(data), phenoData=pd)
source("/nfs/users/nfs_t/ta6/R-Scripts/Ensembl_Stuff.R")
fData(kolo)$feature_symbol <- map_symbol_ensg(rownames(data), is.org="Mmus", is.name="ensg")

saveRDS(kolo, "kolo.rds")

