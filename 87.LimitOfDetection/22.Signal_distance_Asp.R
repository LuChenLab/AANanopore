setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(IRanges)

meta <- as.data.table(openxlsx::read.xlsx("data/MetaInfomation/AtandardAA_LOD.xlsx", sheet = 1))
meta <- meta[amino_acid == "Asp"]
meta[, sig_file := paste0("./analysis/87.LimitOfDetection/01.SelectedL0/", file_id, ".MainL0.txt")]
meta <- meta[file.exists(sig_file)]

bgfile <- meta[concentration == 0, sig_file]
bgsignals <- do.call(rbind, lapply(bgfile, fread))
bgsignals[, ID := paste0("Noise_", ID)]
bgsignals <- bgsignals[Blockade < 0.3 & DwellTime > 0.35]

B_file <- unique(meta[concentration == 0, file_name])
B_file <- paste0("./analysis/81.ABFProcessing/FeatureMatrix/FeatureMatrix_", B_file, ".txt")
Big_FeatureMatrix <- do.call(rbind, lapply(B_file, fread))
Big_FeatureMatrix[, ID := paste0("Noise_", ID)]
Big_FeatureMatrix <- Big_FeatureMatrix[ID %in% bgsignals$ID]

bgsignals0 <- do.call(rbind, lapply(list.files("./analysis/81.ABFProcessing/SelectedSignals", "_background.txt", full.names = TRUE), fread))[Blockade < 0.3 & DwellTime > 0.3]
bgsignals0[, ID := paste0("Noise_", ID)]
Big_FeatureMatrixs <- readRDS("./analysis/81.ABFProcessing/SelectedSignals/01.StandardAA/01.SignalDistance/StandardAA_Big_FeatureMatrixs.Rds")
Big_FeatureMatrixs <- Big_FeatureMatrixs[ID %in% bgsignals0$ID]
Big_FeatureMatrixs[, AA := NULL]

sig0 <- lapply(meta[concentration > 0, sig_file], fread)
sig0 <- do.call(rbind, sig0)[Blockade < 0.3 & DwellTime > 0.35]
B_file <- unique(meta[concentration > 0, file_name])
B_file <- paste0("./analysis/81.ABFProcessing/FeatureMatrix/FeatureMatrix_", B_file, ".txt")
Sig_FeatureMatrix <- do.call(rbind, lapply(B_file, fread))
Sig_FeatureMatrix <- Sig_FeatureMatrix[ID %in% sig0$ID]

if(nrow(Sig_FeatureMatrix) > nrow(Big_FeatureMatrix)) {
  Big_FeatureMatrix_tu <- rbind(Big_FeatureMatrix, Big_FeatureMatrixs[sample(.N, nrow(Sig_FeatureMatrix) - nrow(Big_FeatureMatrix))])
} else {
  Big_FeatureMatrix_tu <- Big_FeatureMatrix[sample(.N, nrow(Sig_FeatureMatrix))]
}

FeatureMatrix <- rbind(Big_FeatureMatrix_tu, Sig_FeatureMatrix)
setkey(FeatureMatrix, ID)

FeatureMatrix <- data.frame(FeatureMatrix[, grepl("^X", colnames(FeatureMatrix)) | colnames(FeatureMatrix) %in% c("DwellTime", "DeltaMean", "StageSD", "CurrentSD", "SignalCurrentPercent", "SignalCurrentWidth", "Blockade"), with = F], row.names = FeatureMatrix[[1]])

dist = as.dist(1 - cor(t(FeatureMatrix), method = "spearman"))
dist <- dbscan::kNN(dist, k = 20)
saveRDS(dist, file = "./analysis/87.LimitOfDetection/02.SignalsDistance/Asp_spearman_distance_knn.Rds")

dist = stats::dist(FeatureMatrix, method = "euclidean")
dist <- dbscan::kNN(dist, k = 20)
saveRDS(dist, file = "./analysis/87.LimitOfDetection/02.SignalsDistance/Asp_euclidean_distance_knn.Rds")







