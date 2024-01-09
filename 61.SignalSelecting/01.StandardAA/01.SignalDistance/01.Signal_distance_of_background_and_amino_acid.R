setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(factoextra)
library(data.table)
library(Biostrings)
library(patchwork)
library(parallel)
library(ggplot2)
library(lsa)

meta1 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 1, cols = 1:7))
meta2 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 2, cols = 1:7))
meta3 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 3, cols = 1:7))
meta4 <- data.table(openxlsx::read.xlsx("./data/单氨基酸实验信息表20210517.xlsx", sheet = 4, cols = 1:7))
meta <- rbind(meta1, meta2, meta3, meta4, use.names = FALSE)
colnames(meta) <- c("file_name", "date", "amino_acid", "concentration", "start_time", "end_time", "type")
meta[amino_acid == "cys", amino_acid := "Cys"]
meta <- meta[concentration == 0 | amino_acid == "blank"]
meta$file_path <- mapply(meta$file_name, FUN = function(x) list.files("./data", recursive = TRUE, full.names = TRUE, pattern = as.character(x))[1])
meta <- meta[!is.na(file_path)]
meta[, start_time := as.numeric(start_time)]
meta[, end_time := as.numeric(end_time)]
# meta <- meta[!grepl("abf", file_name)]
meta[, file_name := gsub(".abf$", "", file_name)]
meta <- meta[file_name %in% gsub("RawSignal_", "", gsub(".MainL0.txt", "", list.files("./analysis/22.SignalSelecting/03.DifferentStates/Background/")))]







# Version 1 ----





mclapply(AMINO_ACID_CODE[1:20], FUN = function(aat){
  print(aat)
  bgsignal <- unique(paste0("./analysis/22.SignalSelecting/03.DifferentStates/Background/RawSignal_", meta[amino_acid == aat, file_name], ".MainL0.txt"))
  bgsignal <- do.call(rbind, lapply(bgsignal, fread))
  bgsignal$File <- stringr::str_remove_all(bgsignal$ID, "_([[:digit:]]+)$")
  
  dataset <- do.call(rbind, lapply(list.files(paste0("./analysis/22.SignalSelecting/03.DifferentStates/", aat), "MainL0.txt", full.names = T), fread))
  dataset$File <- stringr::str_remove_all(dataset$ID, "_([[:digit:]]+)$")
  
  Sig_FeatureMatrix <- paste0("./analysis/21.ABFProcessing/01.StandardAA/FeatureMatrix/FeatureMatrix_", unique(dataset$File), ".txt")
  Big_FeatureMatrix <- paste0("./analysis/21.ABFProcessing/02.Background/FeatureMatrix/FeatureMatrix_", unique(bgsignal$File), ".txt")
  
  Sig_FeatureMatrix <- do.call(rbind, lapply(Sig_FeatureMatrix, fread))
  Big_FeatureMatrix <- do.call(rbind, lapply(Big_FeatureMatrix, fread))
  
  Sig_FeatureMatrix[, ID := paste0(aat, "_", ID)]
  Big_FeatureMatrix[, ID := paste0("Noise_", ID)]
  
  dataset[, ID := paste0(aat, "_", ID)]
  bgsignal[, ID := paste0("Noise_", ID)]
  
  FeatureMatrix <- rbind(Big_FeatureMatrix, Sig_FeatureMatrix)
  setkey(FeatureMatrix, ID)
  FeatureMatrix <- FeatureMatrix[c(dataset$ID, bgsignal$ID), ]
  FeatureMatrix <- data.frame(FeatureMatrix[, grepl("^X", colnames(FeatureMatrix)) | grepl("DwellTime", colnames(FeatureMatrix)), with = F], row.names = FeatureMatrix[[1]])
  
  euclid.dist <- stats::dist(FeatureMatrix, method = "euclidean")
  cosine.dist <- as.dist(1 - lsa::cosine(as.matrix(t(FeatureMatrix))))
  
  res <- list(euclidean = stats::dist(FeatureMatrix, method = "euclidean"), 
              maximum = stats::dist(FeatureMatrix, method = "maximum"), 
              manhattan = stats::dist(FeatureMatrix, method = "manhattan"), 
              canberra = stats::dist(FeatureMatrix, method = "canberra"), 
              binary = stats::dist(FeatureMatrix, method = "binary"), 
              minkowski = stats::dist(FeatureMatrix, method = "minkowski"), 
              cosine = as.dist(1 - lsa::cosine(as.matrix(t(FeatureMatrix)))))
  
  outfile <- file.path("analysis/61.SignalSelecting/01.StandardAA/01.SignalDistance", paste0(aat, ".signal.dist.Rds"))
  saveRDS(res, file = outfile)
}, mc.cores = 10)


aat <- "CbC"

bgsignal <- list.files("./analysis/22.SignalSelecting/03.DifferentStates/CbC", "Background_RawSignal", full.names = T)
bgsignal <- do.call(rbind, lapply(bgsignal, fread))
bgsignal$File <- stringr::str_remove_all(bgsignal$ID, "_([[:digit:]]+)$")

dataset <- do.call(rbind, lapply(list.files(paste0("./analysis/22.SignalSelecting/03.DifferentStates/", aat), "CbC_RawSignal", full.names = T), fread))
dataset$File <- stringr::str_remove_all(dataset$ID, "_([[:digit:]]+)$")

Sig_FeatureMatrix <- paste0("./analysis/21.ABFProcessing/03.MixtureAA/FeatureMatrix/FeatureMatrix_", unique(dataset$File), ".txt")
Big_FeatureMatrix <- paste0("./analysis/21.ABFProcessing/05.CarboxymethylCys/FeatureMatrix/FeatureMatrix_", unique(bgsignal$File), ".txt")

Sig_FeatureMatrix <- do.call(rbind, lapply(Sig_FeatureMatrix, fread))
Big_FeatureMatrix <- do.call(rbind, lapply(Big_FeatureMatrix, fread))

Sig_FeatureMatrix[, ID := paste0(aat, "_", ID)]
Big_FeatureMatrix[, ID := paste0("Noise_", ID)]

dataset[, ID := paste0(aat, "_", ID)]
bgsignal[, ID := paste0("Noise_", ID)]

FeatureMatrix <- rbind(Big_FeatureMatrix, Sig_FeatureMatrix)
setkey(FeatureMatrix, ID)
FeatureMatrix <- FeatureMatrix[c(dataset$ID, bgsignal$ID), ]
FeatureMatrix <- data.frame(FeatureMatrix[, grepl("^X", colnames(FeatureMatrix)) | grepl("DwellTime", colnames(FeatureMatrix)), with = F], row.names = FeatureMatrix[[1]])

res <- list(euclidean = stats::dist(FeatureMatrix, method = "euclidean"), 
            maximum = stats::dist(FeatureMatrix, method = "maximum"), 
            manhattan = stats::dist(FeatureMatrix, method = "manhattan"), 
            canberra = stats::dist(FeatureMatrix, method = "canberra"), 
            binary = stats::dist(FeatureMatrix, method = "binary"), 
            minkowski = stats::dist(FeatureMatrix, method = "minkowski"), 
            cosine = as.dist(1 - lsa::cosine(as.matrix(t(FeatureMatrix)))))

outfile <- file.path("analysis/61.SignalSelecting/01.StandardAA/01.SignalDistance", paste0(aat, ".signal.dist.Rds"))
saveRDS(res, file = outfile)









# Version 2 ----







bgsignals <- mclapply(c(AMINO_ACID_CODE[1:20], "CbC"), function(aat) {
  if(aat != "CbC") {
    bgsignal <- unique(paste0("./analysis/22.SignalSelecting/03.DifferentStates/Background/RawSignal_", meta[amino_acid == aat, file_name], ".MainL0.txt"))
  } else {
    bgsignal <- list.files("./analysis/22.SignalSelecting/03.DifferentStates/CbC", "Background_RawSignal", full.names = T)
  }
  bgsignal <- do.call(rbind, lapply(bgsignal, fread))
  bgsignal$File <- stringr::str_remove_all(bgsignal$ID, "_([[:digit:]]+)$")
  bgsignal[, ID := paste0("Noise_", ID)]
  data.table(bgsignal, AA = aat)
})
bgsignals <- do.call(rbind, bgsignals)



B_file <- bgsignals[, .(File = list(unique(File))), AA]
Big_FeatureMatrixs <- lapply(seq_len(nrow(B_file)), FUN = function(i) {
  Big_FeatureMatrix <- paste0("./analysis/21.ABFProcessing/02.Background/FeatureMatrix/FeatureMatrix_", B_file[[2]][[i]], ".txt")
  if (!all(file.exists(Big_FeatureMatrix))) {
    Big_FeatureMatrix <- paste0("./analysis/21.ABFProcessing/05.CarboxymethylCys/FeatureMatrix/FeatureMatrix_", B_file[[2]][[i]], ".txt")
  }
  Big_FeatureMatrix <- do.call(rbind, lapply(Big_FeatureMatrix, fread))
  Big_FeatureMatrix[, AA := B_file[[1]][i]]
  Big_FeatureMatrix[, ID := paste0("Noise_", ID)]
  Big_FeatureMatrix[ID %in% bgsignals$ID]
})
Big_FeatureMatrixs <- do.call(rbind, Big_FeatureMatrixs)
saveRDS(Big_FeatureMatrixs, "./analysis/61.SignalSelecting/01.StandardAA/01.SignalDistance/StandardAA_Big_FeatureMatrixs.Rds")

mclapply(c(AMINO_ACID_CODE[1:20], "CbC"), FUN = function(aat){
  outfile <- file.path("analysis/61.SignalSelecting/01.StandardAA/01.SignalDistance", paste0(aat, ".signal.dist_V2.Rds"))
  if(file.exists(outfile)) return(NULL)
  print(aat)
  if(aat == "CbC") {
    dataset <- do.call(rbind, lapply(list.files(paste0("./analysis/22.SignalSelecting/03.DifferentStates/", aat), "CbC_RawSignal", full.names = T), fread))
    dataset$File <- stringr::str_remove_all(dataset$ID, "_([[:digit:]]+)$")
    Sig_FeatureMatrix <- paste0("./analysis/21.ABFProcessing/03.MixtureAA/FeatureMatrix/FeatureMatrix_", unique(dataset$File), ".txt")
  } else {
    dataset <- do.call(rbind, lapply(list.files(paste0("./analysis/22.SignalSelecting/03.DifferentStates/", aat), "MainL0.txt", full.names = T), fread))
    dataset$File <- stringr::str_remove_all(dataset$ID, "_([[:digit:]]+)$")
    Sig_FeatureMatrix <- paste0("./analysis/21.ABFProcessing/01.StandardAA/FeatureMatrix/FeatureMatrix_", unique(dataset$File), ".txt")
  }
  Sig_FeatureMatrix <- do.call(rbind, lapply(Sig_FeatureMatrix, fread))
  Sig_FeatureMatrix <- Sig_FeatureMatrix[ID %in% dataset[, ID]]
  Sig_FeatureMatrix[, ID := paste0(aat, "_", ID)]
  Sig_FeatureMatrix[, AA := aat]
  
  if(Big_FeatureMatrixs[AA == aat, .N] >= Sig_FeatureMatrix[, .N]) {
    Big_FeatureMatrix <- rbind(Big_FeatureMatrixs[AA == aat, .SD[sample(.N, Sig_FeatureMatrix[, .N])], ], 
                               Big_FeatureMatrixs[AA != aat, .SD[sample(.N, Sig_FeatureMatrix[, .N])], ])
    
  } else {
    Big_FeatureMatrix <- rbind(Big_FeatureMatrixs[AA == aat, .SD], 
                               Big_FeatureMatrixs[AA != aat, .SD[sample(.N, Sig_FeatureMatrix[, .N] * 2 - Big_FeatureMatrixs[AA == aat, .N]), ]])
  }
  
  FeatureMatrix <- rbind(Big_FeatureMatrix, Sig_FeatureMatrix)
  setkey(FeatureMatrix, ID)
  FeatureMatrix <- data.frame(FeatureMatrix[, grepl("^X", colnames(FeatureMatrix)) | grepl("DwellTime", colnames(FeatureMatrix)), with = F], row.names = FeatureMatrix[[1]])
  
  res <- list(euclidean = stats::dist(FeatureMatrix, method = "euclidean"),
              maximum = stats::dist(FeatureMatrix, method = "maximum"),
              manhattan = stats::dist(FeatureMatrix, method = "manhattan"),
              canberra = stats::dist(FeatureMatrix, method = "canberra"),
              binary = stats::dist(FeatureMatrix, method = "binary"),
              minkowski = stats::dist(FeatureMatrix, method = "minkowski"),
              cosine = as.dist(1 - lsa::cosine(as.matrix(t(FeatureMatrix)))))

  saveRDS(res, file = outfile)
  dim(FeatureMatrix)
}, mc.cores = 11)
