setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(factoextra)
library(data.table)
library(Biostrings)
library(patchwork)
library(parallel)
library(ggplot2)
library(lsa)


bgfile <- list.files("./analysis/82.ABFProcessing/SelectedSignals", "_background.txt", full.names = TRUE)
bgsignals <- do.call(rbind, lapply(bgfile, fread))
bgsignals[, ID := paste0("Noise_", ID)]

B_file <- bgsignals[, .(File = list(unique(file_id))), AA]
Big_FeatureMatrixs <- lapply(seq_len(nrow(B_file)), FUN = function(i) {
  FMfs <- mapply(function(x) x[1], strsplit(B_file[[2]][[i]], "\\."))
  Big_FeatureMatrix <- paste0("./analysis/82.ABFProcessing/FeatureMatrix/FeatureMatrix_", FMfs, ".txt")
  if (!all(file.exists(Big_FeatureMatrix))) return(NULL)
  Big_FeatureMatrix <- do.call(rbind, lapply(Big_FeatureMatrix, fread))
  Big_FeatureMatrix[, AA := B_file[[1]][i]]
  Big_FeatureMatrix[, ID := paste0("Noise_", ID)]
  Big_FeatureMatrix[ID %in% bgsignals$ID]
})
Big_FeatureMatrixs <- do.call(rbind, Big_FeatureMatrixs)
saveRDS(Big_FeatureMatrixs, "./analysis/82.ABFProcessing/SelectedSignals/01.StandardAA/01.SignalDistance/StandardAA_Big_FeatureMatrixs.Rds")



mclapply(c(AMINO_ACID_CODE[1:20], "CbC"), FUN = function(aat){
  print(aat)
  outfile <- file.path("./analysis/82.ABFProcessing/SelectedSignals/01.StandardAA/01.SignalDistance", paste0(aat, ".signal.dist.Rds"))
  if(file.exists(outfile)) return(NULL)
  
  dataset <- fread(paste0("./analysis/82.ABFProcessing/SelectedSignals/", aat, "_signal.txt"))
  dataset[, file_id := as.character(file_id)]
  FMfs <- mapply(function(x) x[1], strsplit(unique(dataset$file_id), "\\."))
  
  Sig_FeatureMatrix <- paste0("./analysis/82.ABFProcessing/FeatureMatrix/FeatureMatrix_", FMfs, ".txt")
  if (!all(file.exists(Sig_FeatureMatrix))) return(NULL)
  
  Sig_FeatureMatrix <- do.call(rbind, lapply(Sig_FeatureMatrix, fread))
  Sig_FeatureMatrix <- Sig_FeatureMatrix[ID %in% dataset[, ID]]
  Sig_FeatureMatrix[, ID := paste0(aat, "_", ID)]
  Sig_FeatureMatrix[, AA := aat]
  
  if(Big_FeatureMatrixs[AA == aat, .N] >= Sig_FeatureMatrix[, .N]) {
    Big_FeatureMatrix <- rbind(Big_FeatureMatrixs[AA == aat, .SD[sample(.N, Sig_FeatureMatrix[, .N])], ], 
                               Big_FeatureMatrixs[AA != aat, .SD[sample(.N, Sig_FeatureMatrix[, .N])], ])
    
  } else {
    Big_FeatureMatrix <- rbind(Big_FeatureMatrixs[AA == aat, .SD], 
                               Big_FeatureMatrixs[AA != aat, .SD[sample(.N, Sig_FeatureMatrix[, .N] * 1 - Big_FeatureMatrixs[AA == aat, .N]), ]])
  }
  
  FeatureMatrix <- rbind(Big_FeatureMatrix, Sig_FeatureMatrix)
  setkey(FeatureMatrix, ID)
  FeatureMatrix <- data.frame(FeatureMatrix[, grepl("^X", colnames(FeatureMatrix)) | grepl("DwellTime", colnames(FeatureMatrix)), with = F], row.names = FeatureMatrix[[1]])
  euclidean = stats::dist(FeatureMatrix, method = "euclidean")
  saveRDS(euclidean, file = outfile)
}, mc.cores = 21)


