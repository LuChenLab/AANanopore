setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(parallel)

files <- list.files("./analysis/11.SignalIdentification/Jan07/AASignal", pattern = "txt", full.names = TRUE)
files <- grep(files, pattern = "Selected", invert = T, value = T)
files <- grep(files, pattern = "Feature_Mat", invert = T, value = T)
# List <- lapply(files, fread)
# names(List) <- gsub("RawSignal_", "", gsub(".txt", "", basename(files)))

Mat_List <- mclapply(seq_along(files), function(i) {
  print(i)
  if(file.exists(paste0("./analysis/13.MachineLearning/01.DataPreparation/FeatureMatrix/Density200Points/StandardAA/", gsub(".txt", "_Feature_Mat.txt", basename(files[i]))))) return(NULL)
  Mat <- fread(files[i])
  Current <- readRDS(paste0("./analysis/11.SignalIdentification/Dec27/SignalCurrent_", gsub("RawSignal_", "", gsub(".txt", "", basename(files[i]))), ".Rds"))
  Feature_Mat <- lapply(Mat$ID, function(j) {
    D <- density(Current[ID == j, Current], from = 0, to = 1, n = 200, adjust = 0.5)$y
    pA <- sort(Current[ID == j, Current])
    pA <- pA[round(seq(1, length(pA), length.out = 100))]
    D <- c(D, pA)
    cbind(t(data.frame(round(D, 4), row.names = c(paste0("X", sprintf("%03d", 1:200)), paste0("P", sprintf("%03d", 1:100))))), 
          Mat[ID == j, .(DeltaMean, StageSD, CurrentSD, Segments, SignalCurrentPercent, L2Ratio, DwellTime, Blockade, A, ID)])
  })
  Feature_Mat <- do.call(rbind, Feature_Mat)
  fwrite(Feature_Mat, paste0("./analysis/13.MachineLearning/01.DataPreparation/FeatureMatrix/Density200Points/StandardAA/", gsub(".txt", "_Feature_Mat.txt", basename(files[i]))), sep = "\t", quote = F)
}, mc.cores = 1)


files2 <- list.files("./analysis/11.SignalIdentification/Jan07/BackgroundSignal", pattern = ".txt", full.names = TRUE)
files2 <- grep(files2, pattern = "Selected", invert = T, value = T)
files2 <- grep(files2, pattern = "Feature_Mat", invert = T, value = T)

Mat_List <- mclapply(seq_along(files2), function(i) {
  print(i)
  Current <- readRDS(paste0("./analysis/11.SignalIdentification/Dec27/SignalCurrent_", gsub("RawSignal_", "", gsub(".txt", "", basename(files2[i]))), ".Rds"))
  Mat <- fread(files2[i])
  Feature_Mat <- lapply(Mat$ID, function(j) {
    D <- density(Current[ID == j, Current], from = 0, to = 1, n = 200, adjust = 0.5)$y
    pA <- sort(Current[ID == j, Current])
    pA <- pA[round(seq(1, length(pA), length.out = 100))]
    D <- c(D, pA)
    cbind(t(data.frame(round(D, 4), row.names = c(paste0("X", sprintf("%03d", 1:200)), paste0("P", sprintf("%03d", 1:100))))), 
          Mat[ID == j, .(DeltaMean, StageSD, CurrentSD, Segments, SignalCurrentPercent, L2Ratio, DwellTime, Blockade, A, ID)])
  })
  Feature_Mat <- do.call(rbind, Feature_Mat)
  if(!is.null(Feature_Mat)) {
    fwrite(Feature_Mat, paste0("./analysis/13.MachineLearning/01.DataPreparation/FeatureMatrix/Density200Points/BackgroundSignal/", gsub(".txt", "_Feature_Mat.txt", basename(files2[i]))), sep = "\t", quote = F)
  }
}, mc.cores = 1)
