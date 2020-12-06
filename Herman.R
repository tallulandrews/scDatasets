# Herman, Sagar, Grun. FateID infers cell fate bias in multipotent progenitors from single-cell RNA-seq data. 
# PMID: 29630061
# GSE100037

mat_files <- c("GSM2668205_LSK1_merged.coutt.csv.gz",
"GSM2668206_LSK2_merged.coutt.csv.gz",
"GSM2668207_LMPP1_1_merged.coutt.csv.gz",
"GSM2668208_LMPP1_2_merged.coutt.csv.gz",
"GSM2668209_LMPP2_1_merged.coutt.csv.gz",
"GSM2668210_LMPP2_2_merged.coutt.csv.gz",
"GSM2668211_CLP1_merged.coutt.csv.gz",
"GSM2668212_CLP2_merged.coutt.csv.gz",
"GSM2668213_M1_p1_1.coutt.csv.gz",
"GSM2668214_M1_p1_2.coutt.csv.gz",
"GSM2668215_M1_p1_3.coutt.csv.gz",
"GSM2668216_M1_p1_4.coutt.csv.gz",
"GSM2668217_M1_p2_1.coutt.csv.gz",
"GSM2668218_M1_p2_2.coutt.csv.gz",
"GSM2668219_M1_p2_34.coutt.csv.gz",
"GSM2668220_M3_p1_1.coutt.csv.gz",
"GSM2668221_M3_p1_2.coutt.csv.gz",
"GSM2668222_M3_p1_3.coutt.csv.gz",
"GSM2668223_M3_p1_4.coutt.csv.gz",
"GSM2668224_NP1_9.coutt.csv.gz",
"GSM2668225_NP1_10.coutt.csv.gz",
"GSM2668226_ES_2_1fold.coutt.csv.gz",
"GSM2668227_ES_2_5fold.coutt.csv.gz",
"GSM2668228_ES_2_7fold.coutt.csv.gz",
"GSM2668229_ES_2_10fold.coutt.csv.gz",
"GSM2857521_CLP_WT3615_IL7rpos_gate13_1.coutt.csv.gz",
"GSM2857522_CLP_WT3615_Cd34pos_gate15_2.coutt.csv.gz",
"GSM2857523_CLP_WT3616_IL7rpos_gate13_3.coutt.csv.gz",
"GSM2857524_CLP_WT3616_Cd34pos_gate15_4.coutt.csv.gz",
"GSM2857525_WT_3615_3616_mixed_1.coutt.csv.gz",
"GSM2857526_WT_3615_3616_mixed_2.coutt.csv.gz",
"GSM2857527_LMPP_WT3615_gate7_3.coutt.csv.gz",
"GSM2857528_LMPP_WT3616_gate7_4.coutt.csv.gz",
"GSM2857529_WT_gate16_cult7d_Flt3L_1.coutt.csv.gz",
"GSM2857530_WT_gate16_cult7d_Flt3L_2.coutt.csv.gz",
"GSM2857531_WT_gate18_cult7d_Flt3L_1.coutt.csv.gz",
"GSM2857532_WT_gate18_cult7d_Flt3L_2.coutt.csv.gz",
"GSM2857533_WT_gate17_cult7d_Flt3L_3.coutt.csv.gz",
"GSM2857534_WT_gate17_cult7d_Flt3L_4.coutt.csv.gz")


out_obj <- vector();
exp1_cells <- read.table("GSE100037_cellbarcodes_cellid.csv.gz", sep=",", header=T, stringsAsFactors=FALSE)
exp2_cells <- read.table("GSE100037_cellbarcodes_cellid_culture_experiment_2.csv.gz", sep=",", header=T, stringsAsFactors=FALSE)

for (f in mat_files) {
	require(Matrix)
	require(SingleCellExperiment)
	thing <- read.table(f, stringsAsFactors=FALSE, header=T)
	rownames(thing) <- thing[,1]
	thing <- thing[,-1]
	thing <- Matrix(as.matrix(thing))
	type <- sub("GSM[1234567890]*_", "", f)
	type <- sub(".coutt.csv.gz", "", type)
	type <- sub("_merged", "", type)
	colnames(thing) <- paste(type, 1:ncol(thing), sep="_")
	exp <- rep(0, length=ncol(thing));
	exp[colnames(thing) %in% exp1_cells$cellid] <- 1;
	exp[colnames(thing) %in% exp2_cells$cellid] <- 2;

	obj <- SingleCellExperiment(assays=list(counts=round(thing)), colData=data.frame(type=rep(type, ncol(thing)), experiment=exp, stringsAsFactors=FALSE))
	if (is.null(dim(out_obj))) {
		out_obj <- obj
	} else {
		merge1 <- as.matrix(out_obj@assays[["counts"]])
		merge2 <- as.matrix(obj@assays[["counts"]])
		genes <- sort(union(rownames(merge1), rownames(merge2)))
		merge1 <- merge1[match(genes, rownames(merge1)),]
		merge2 <- merge2[match(genes, rownames(merge2)),]
		rownames(merge1) <- genes; rownames(merge2) <- genes;
		merge1[is.na(merge1)] <- 0
		merge2[is.na(merge2)] <- 0
		merge1 <- Matrix(merge1)
		merge2 <- Matrix(merge2)
		out_obj <- SingleCellExperiment(assays=list(counts=cbind(merge1, merge2)), colData=rbind(colData(out_obj), colData(obj)))
	}
}

saveRDS(out_obj, "Herman.rds")

#out.qc <- out_obj[, Matrix::colSums(out_obj@assays[["counts"]]) >= 2000]
#qc2 <- out.qc@assays[["counts"]][rownames(out.qc) == "Kcnq1ot1",]/Matrix::colSums(out.qc@assays[["counts"]]);
#out.qc <- out.qc[,qc2 <= 0.02]
#out.qc <- out.qc[Matrix::rowSums(out.qc@assays[["counts"]]) > 100,]
#
#set.seed(4201)
#library(RaceID)
#sc <- SCseq(out.qc@assays[["counts"]])
#sc <- filterdata(sc,mintotal=2000)
#sc <- compdist(sc,metric="pearson")
#sc <- clustexp(sc)
#sc <- findoutliers(sc)
#sc <- compfr(sc,knn=50)


#mintotal = 2000, minexpr = 3, outminc = 3, FSelect = TRUE, probthr = 10−4 and random-forest-based reclassification
