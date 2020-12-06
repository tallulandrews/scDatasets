# array data for hiPSC-derived cholangiocyte-like cells 
# Sampaziotis et al. (2015). Cholangiocytes derived from human induced pluripotent stem cells for disease modeling and drug validation. Nature Biotechnology. 33:845-852
# PMID: 26167629
# ArrayExpress: E-MTAB-2965
#require("lumi") # I think this is already QCed & quantile normalized

probe_lvls <- read.delim("fotis_sampaziotis_2014-09-15_groutset_quantile_nobkgd_ArrayExpress_DataFile.txt", header=TRUE, stringsAsFactors=FALSE)

probe_lvls <- probe_lvls[-c(1:6),]
colnames(probe_lvls) <- probe_lvls[1,]
probe_lvls <- probe_lvls[-1,]

sample_anno <- read.delim("E-MTAB-2965.sdrf.txt", header=TRUE, stringsAsFactors=FALSE)
probe_anno <- read.delim("A-GEOD-10558.adf.txt", header=TRUE, stringsAsFactors=FALSE)


library("illuminaHumanv4.db")
thing <- mget(x=probe_lvls[,1], envir= illuminaHumanv4SYMBOL)
probe_lvls$Gene <- unlist(thing)

require("limma")
limma_obj <- read.ilmn("fotis_sampaziotis_2014-09-15_groutset_quantile_nobkgd_ArrayExpress_DataFile.txt")
conditions <- sample_anno[match(colnames(limma_obj), sample_anno[,1]), c("Source.Name", "Description")]
conditions <- matrix(unlist(strsplit(conditions[,2], " ")), ncol=2,byrow=TRUE)

probe_filter <- rowSums(limma_obj$other[[1]] < 0.05) >= 3
limma_obj<- limma_obj[probe_filter,] 
cond <- factor(conditions[,1], levels=c("hPSCs", "HBs", "CLCs", "BIL"))
mod <- model.matrix(~cond)

fit <- lmFit(limma_obj$E, mod)
contrast.matrix <- makeContrasts(condCLCs-condBIL, levels=mod)
fit_immature_vs_mature_Chol <- contrasts.fit(fit, contrast.matrix)
fit_immature_vs_mature_Chol <- eBayes(fit_immature_vs_mature_Chol)
results_imm_mat <- decideTests(fit_immature_vs_mature_Chol)

# List of genes I want: 
#	CLCs & BIL vs hPSCs & HBs = Chol
#	CLCs vs hPSCs & HBs OR BIL vs hPSCs & HBs = Chol
#	CLCs vs BIL & in one of the above = immature/mature Chol

contrast.matrix <- makeContrasts(condCLCs-condBIL, condCLCs+condBIL-Intercept+condHBs, condCLCs-Intercept, condCLCs-condHBs, condBIL-Intercept, condBIL-condHBs, levels=mod)
fit_contrasts <- contrasts.fit(fit, contrast.matrix)
fit_contrasts <- eBayes(fit_contrasts)

results <- decideTests(fit_contrasts)
Chol_markers <- (results[,3]==1 & results[,4]==1) | (results[,5]==1 & results[,6]==1)
young_chol_markers <- results_imm_mat > 0
mature_chol_markers <- results_imm_mat < 0

Mature_score <- abs(fit_contrasts$coefficients[mature_chol_markers & Chol_markers,1] ) + fit_contrasts$coefficients[mature_chol_markers & Chol_markers,5] + fit_contrasts$coefficients[mature_chol_markers & Chol_markers,6]
Mature_score <- Mature_score/3

Immature_score <- fit_contrasts$coefficients[young_chol_markers & Chol_markers,1] + fit_contrasts$coefficients[young_chol_markers & Chol_markers,3] + fit_contrasts$coefficients[young_chol_markers & Chol_markers,4]
Immature_score <- Immature_score/3

AllChol_score <- rowMeans(fit_contrasts$coefficients[Chol_markers & !young_chol_markers & !mature_chol_markers,c(3:6)])

Mature_score<-Mature_score[order(-Mature_score)]
Immature_score<-Immature_score[order(-Immature_score)]
AllChol_score<-AllChol_score[order(-AllChol_score)]

library("illuminaHumanv4.db")
symbol <-unlist(mget(x=names(Mature_score), envir=illuminaHumanv4SYMBOL))
Mature_genes <- Mature_score[!is.na(symbol)]
names(Mature_genes) <- symbol[!is.na(symbol)]

symbol <-unlist(mget(x=names(Immature_score), envir=illuminaHumanv4SYMBOL))
Immature_genes <- Immature_score[!is.na(symbol)]
names(Immature_genes) <- symbol[!is.na(symbol)]

symbol <-unlist(mget(x=names(AllChol_score), envir=illuminaHumanv4SYMBOL))
AllChol_genes <- AllChol_score[!is.na(symbol)]
names(AllChol_genes) <- symbol[!is.na(symbol)]

OUT <- list(Chol=AllChol_genes, CLC=Immature_genes, Mature=Mature_genes)

saveRDS(OUT, file="Sampaziotis_CholMarkers.rds")
