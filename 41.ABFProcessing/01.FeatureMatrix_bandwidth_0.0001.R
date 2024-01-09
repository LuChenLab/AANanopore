setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")

library(data.table)
library(Biostrings)
library(IRanges)
library(ggplot2)
library(parallel)
library(changepoint)

Sig_files <- list.files('./analysis/21.ABFProcessing/01.StandardAA/RawSignal', 'RawSignal_', full.names = T)
Cur_files <- list.files('./analysis/21.ABFProcessing/01.StandardAA/SignalCurrent', 'SignalCurrent_', full.names = T)
stopifnot(identical(gsub('.txt', '', gsub('RawSignal_', '', basename(Sig_files))), gsub('.Rds', '', gsub('SignalCurrent_', '', basename(Cur_files)))))

mclapply(seq_along(Sig_files), function(i) {
  print(i)
  Sigs <- fread(Sig_files[i])
  Current <- readRDS(Cur_files[i])
  
  # Feature matrix
  Mat <- lapply(Sigs[, ID], function(j) {
    D <- density(Current[ID == j, Current], from = 0, to = 1, n = 1024, bw = 1e-04)$y
    pA <- Current[ID == j, Current]
    pA <- pA[round(seq(1, length(pA), length.out = 512))]
    c(D, pA)
  })
  Mat <- as.data.table(do.call(rbind, Mat))
  colnames(Mat) <- c(paste0("X", sprintf("%04d", seq_len(1024))), paste0("P", sprintf("%03d", 1:512)))
  Mat <- data.table(ID = Sigs[, ID], Mat)
  Mat <- merge(Mat, Sigs, by = "ID")
  setkey(Mat, ID)
  Mat <- Mat[Sigs$ID, ]
  fwrite(Mat, paste0("./analysis/41.ABFProcessing/01.StandardAA/FeatureMatrix/FeatureMatrix_", gsub('RawSignal_', '', basename(Sig_files[i]))), sep = "\t", quote = F)
}, mc.cores = 10)





Sig_files <- list.files('./analysis/21.ABFProcessing/02.Background/RawSignal', 'RawSignal_', full.names = T)
Cur_files <- list.files('./analysis/21.ABFProcessing/02.Background/SignalCurrent', 'SignalCurrent_', full.names = T)
stopifnot(identical(gsub('.txt', '', gsub('RawSignal_', '', basename(Sig_files))), gsub('.Rds', '', gsub('SignalCurrent_', '', basename(Cur_files)))))

mclapply(seq_along(Sig_files), function(i) {
  print(i)
  Sigs <- fread(Sig_files[i])
  Current <- readRDS(Cur_files[i])
  
  # Feature matrix
  Mat <- lapply(Sigs[, ID], function(j) {
    D <- density(Current[ID == j, Current], from = 0, to = 1, n = 1024, bw = 1e-04)$y
    pA <- Current[ID == j, Current]
    pA <- pA[round(seq(1, length(pA), length.out = 512))]
    c(D, pA)
  })
  Mat <- as.data.table(do.call(rbind, Mat))
  colnames(Mat) <- c(paste0("X", sprintf("%04d", seq_len(1024))), paste0("P", sprintf("%03d", 1:512)))
  Mat <- data.table(ID = Sigs[, ID], Mat)
  Mat <- merge(Mat, Sigs, by = "ID")
  setkey(Mat, ID)
  Mat <- Mat[Sigs$ID, ]
  fwrite(Mat, paste0("./analysis/41.ABFProcessing/02.Background/FeatureMatrix/FeatureMatrix_", gsub('RawSignal_', '', basename(Sig_files[i]))), sep = "\t", quote = F)
}, mc.cores = 10)






Sig_files <- list.files('./analysis/21.ABFProcessing/03.MixtureAA/RawSignal', 'RawSignal_', full.names = T)
Cur_files <- list.files('./analysis/21.ABFProcessing/03.MixtureAA/SignalCurrent', 'SignalCurrent_', full.names = T)
stopifnot(identical(gsub('.txt', '', gsub('RawSignal_', '', basename(Sig_files))), gsub('.Rds', '', gsub('SignalCurrent_', '', basename(Cur_files)))))

lapply(seq_along(Sig_files), function(i) {
  print(i)
  if(file.exists(paste0("./analysis/41.ABFProcessing/03.MixtureAA/FeatureMatrix/FeatureMatrix_", gsub('RawSignal_', '', basename(Sig_files[i]))))) return(NULL)
  Sigs <- fread(Sig_files[i])
  Current <- readRDS(Cur_files[i])
  
  # Feature matrix
  Mat <- lapply(Sigs[, ID], function(j) {
    D <- density(Current[ID == j, Current], from = 0, to = 1, n = 1024, bw = 1e-04)$y
    pA <- Current[ID == j, Current]
    pA <- pA[round(seq(1, length(pA), length.out = 512))]
    c(D, pA)
  })
  Mat <- as.data.table(do.call(rbind, Mat))
  colnames(Mat) <- c(paste0("X", sprintf("%04d", seq_len(1024))), paste0("P", sprintf("%03d", 1:512)))
  Mat <- data.table(ID = Sigs[, ID], Mat)
  Mat <- merge(Mat, Sigs, by = "ID")
  setkey(Mat, ID)
  Mat <- Mat[Sigs$ID, ]
  fwrite(Mat, paste0("./analysis/41.ABFProcessing/03.MixtureAA/FeatureMatrix/FeatureMatrix_", gsub('RawSignal_', '', basename(Sig_files[i]))), sep = "\t", quote = F)
})

