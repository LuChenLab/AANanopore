setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(IRanges)

meta <- fread("./data/ChenShanchuan/20240103/meta_20240103.txt", sep = "\t")
meta[, sig_file := paste0("./analysis/93.Revision/02.PolypeptideHydrolysis/02.SelectedL0/", file_id, ".MainL0.txt")]
meta[amino_acid == "", amino_acid := NA]
meta0 <- meta[file.exists(sig_file)]

bgfile <- meta0[is.na(amino_acid), sig_file]
bgsignals <- do.call(rbind, lapply(bgfile, fread))
bgsignals[, ID := paste0("Noise_", ID)]
bgsignals <- bgsignals[Blockade < 0.3 & DwellTime > 0.3]

B_file <- unique(meta0[is.na(amino_acid), file_name])
B_file <- paste0("./analysis/81.ABFProcessing/FeatureMatrix/FeatureMatrix_", B_file, ".txt")
Big_FeatureMatrix <- do.call(rbind, lapply(B_file, fread))
Big_FeatureMatrix[, ID := paste0("Noise_", ID)]
Big_FeatureMatrix <- Big_FeatureMatrix[ID %in% bgsignals$ID]

bgsignals0 <- do.call(rbind, lapply(list.files("./analysis/81.ABFProcessing/SelectedSignals", "_background.txt", full.names = TRUE), fread))[Blockade < 0.3 & DwellTime > 0.3]
bgsignals0[, ID := paste0("Noise_", ID)]
Big_FeatureMatrixs <- readRDS("./analysis/81.ABFProcessing/SelectedSignals/01.StandardAA/01.SignalDistance/StandardAA_Big_FeatureMatrixs.Rds")
Big_FeatureMatrixs <- Big_FeatureMatrixs[ID %in% bgsignals0$ID]
Big_FeatureMatrixs[, AA := NULL]

sig0 <- lapply(meta0[!is.na(amino_acid), sig_file], fread)
sig0 <- do.call(rbind, sig0)[Blockade < 0.3 & DwellTime > 0.35]
B_file <- unique(meta0[!is.na(amino_acid), file_name])
B_file <- paste0("./analysis/81.ABFProcessing/FeatureMatrix/FeatureMatrix_", B_file, ".txt")
Sig_FeatureMatrix <- do.call(rbind, lapply(B_file, fread))
Sig_FeatureMatrix <- Sig_FeatureMatrix[ID %in% sig0$ID]

# if(nrow(Sig_FeatureMatrix) > nrow(Big_FeatureMatrix)) {
#   Big_FeatureMatrix_tu <- rbind(Big_FeatureMatrix, Big_FeatureMatrixs[sample(.N, nrow(Sig_FeatureMatrix) - nrow(Big_FeatureMatrix))])
# } else {
#   Big_FeatureMatrix_tu <- Big_FeatureMatrix[sample(.N, nrow(Sig_FeatureMatrix))]
# }
# 
# FeatureMatrix <- rbind(Big_FeatureMatrix_tu, Sig_FeatureMatrix)
# setkey(FeatureMatrix, ID)
# mean(FeatureMatrix[, grepl("^Noise", ID)])
# 
# FeatureMatrix <- data.frame(FeatureMatrix[, grepl("^X", colnames(FeatureMatrix)) | colnames(FeatureMatrix) %in% c("DwellTime", "DeltaMean", "StageSD", "CurrentSD", "SignalCurrentPercent", "SignalCurrentWidth", "Blockade"), with = F], row.names = FeatureMatrix[[1]])
# stopifnot(mean(grepl("^Noise", row.names(FeatureMatrix))) == 0.5)
# dim(FeatureMatrix)
# 
# dist = stats::dist(FeatureMatrix, method = "euclidean")
# dist <- dbscan::kNN(dist, k = 20)
# saveRDS(dist, file = "./analysis/93.Revision/02.PolypeptideHydrolysis/03.SignalDistance/euclidean_distance_knn.Rds")
# 

Big_FeatureMatrix <- data.table()
mclapply(sig0[, unique(A)], function(x) {
  print(x)
  outfile <- paste0("./analysis/93.Revision/02.PolypeptideHydrolysis/03.SignalDistance/", x)
  if(file.exists(paste0(outfile, "_spearman_distance_knn.Rds")) & file.exists(paste0(outfile, "_euclidean_distance_knn.Rds"))) return(NULL)
  sig1 <- sig0[A == x]
  Sig_FM <- Sig_FeatureMatrix[ID %in% sig1$ID]
  
  if(nrow(Sig_FM) > nrow(Big_FeatureMatrix)) {
    Big_FeatureMatrix_tu <- rbind(Big_FeatureMatrix, Big_FeatureMatrixs[sample(.N, nrow(Sig_FM) - nrow(Big_FeatureMatrix))])
  } else {
    Big_FeatureMatrix_tu <- Big_FeatureMatrix[sample(.N, nrow(Sig_FM))]
  }
  
  FeatureMatrix <- rbind(Big_FeatureMatrix_tu, Sig_FM)
  setkey(FeatureMatrix, ID)
  stopifnot(mean(FeatureMatrix[, grepl("^Noise", ID)]) == 0.5)
  
  FeatureMatrix <- data.frame(FeatureMatrix[, grepl("^X", colnames(FeatureMatrix)) | colnames(FeatureMatrix) %in% c("DwellTime", "DeltaMean", "StageSD", "CurrentSD", "SignalCurrentPercent", "SignalCurrentWidth", "Blockade"), with = F], row.names = FeatureMatrix[[1]])
  stopifnot(mean(grepl("^Noise", row.names(FeatureMatrix))) == 0.5)
  
  # dist <- as.dist(1 - cor(t(FeatureMatrix), method = "spearman"))
  # dist <- dbscan::kNN(dist, k = 20)
  # saveRDS(dist, file = paste0(outfile, "_spearman_distance_knn.Rds"))
  
  dist = stats::dist(FeatureMatrix, method = "euclidean")
  dist <- dbscan::kNN(dist, k = 20)
  saveRDS(dist, file = paste0(outfile, "_euclidean_distance_knn.Rds"))
}, mc.cores = length(sig0[, unique(A)]))

