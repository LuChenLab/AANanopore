setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")

library(data.table)
library(Biostrings)
library(IRanges)
library(ggplot2)
library(parallel)
library(changepoint)

MainRidge <- function(x, bw = 0.005, n = 1024, plot = F, ...) {
  bw <- bw * max(abs(x))
  den <- density(x, bw = bw, n = n, ...)
  MainPeak <- with(den, x[which.max(y)])
  uR <- IRanges(diff(den$y) >= 0)
  dR <- IRanges(diff(den$y) <= 0)
  while (length(uR) != length(dR)) {
    bw <- bw * 1.05
    den <- density(x, bw = bw, n = n, ...)
    MainPeak <- with(den, x[which.max(y)])
    uR <- IRanges(diff(den$y) >= 0)
    dR <- IRanges(diff(den$y) <= 0)
  }
  res <- do.call(c, lapply(seq_along(uR), function(i) {
    reduce(c(uR[i], dR[i]))
  }))
  res <- res[subjectHits(findOverlaps(IRanges(which.max(den$y), which.max(den$y)), res))]
  MainRidgeRatio <- mean(x > den$x[start(res)] & x < den$x[end(res)])
  while (MainRidgeRatio < 0.1) {
    bw <- bw * 1.05
    den <- density(x, bw = bw, n = n, ...)
    MainPeak <- with(den, x[which.max(y)])
    uR <- IRanges(diff(den$y) >= 0)
    dR <- IRanges(diff(den$y) <= 0)
    while (length(uR) != length(dR)) {
      bw <- bw * 1.05
      den <- density(x, bw = bw, n = n, ...)
      MainPeak <- with(den, x[which.max(y)])
      uR <- IRanges(diff(den$y) >= 0)
      dR <- IRanges(diff(den$y) <= 0)
    }
    res <- do.call(c, lapply(seq_along(uR), function(i) {
      reduce(c(uR[i], dR[i]))
    }))
    res <- res[subjectHits(findOverlaps(IRanges(which.max(den$y), which.max(den$y)), res))]
    MainRidgeRatio <- mean(x > den$x[start(res)] & x < den$x[end(res)])
  }
  MainRidgeWidth <- den$x[end(res)] - den$x[start(res)]
  if(plot) {
    plot(den)
    abline(v = with(den, x[c(start(res), end(res))]))
    abline(v = MainPeak, col = 2)
  }
  data.table(MainPeak, MainRidgeWidth, MainRidgeRatio)
}


KBins2 <- function(sig, minseglen1 = 100, minseglen2 = 10, pen.value = 1e-3) {
  ansmean1 <- suppressWarnings(changepoint::cpt.mean(sig, penalty = "MBIC", method = "PELT", minseglen = minseglen1))
  Tab1 <- data.table::data.table(P = sig, B1 = rep(seq_len(changepoint::nseg(ansmean1)), changepoint::seg.len(ansmean1)))
  Tab2 <- Tab1[, .(P = median(P)), by = "B1"]
  
  ansmean2 <- suppressWarnings(changepoint::cpt.meanvar(Tab2[, P], penalty = "Asymptotic", pen.value = pen.value, method = "PELT", minseglen = minseglen2))
  Tab2[, B := rep(seq_len(changepoint::nseg(ansmean2)), changepoint::seg.len(ansmean2))]
  res <- merge(Tab1, Tab2[, .(B1, B)], by = "B1")
  res[, B1 := NULL]
  res[, .(P = median(P), N = .N), "B"][, rep(P, N)]
}

mean2 <- function(x, q1 = 0.25, q2 = 0.75) {
  q <- quantile(x, c(q1, q2))
  mean(x[x > min(q) & x < max(q)])
}

WhiskerRange <- function(x) {
  iqr <- 1.5 * IQR(x, na.rm = TRUE)
  c(min(x[x >= quantile(x, 1/4, na.rm = T) - iqr], na.rm = TRUE), max(x[x <= quantile(x, 3/4, na.rm = T) + iqr], na.rm = TRUE))
}

GetSignal <- function(x, minblockade = 0.1) {
  abf2 <- x[rep(runLength(Rle(x$Sm)) >= 20, runLength(Rle(x$Sm)))]
  gr <- IRanges(end = cumsum(runLength(Rle(abf2[, Sm]))), width = runLength(Rle(abf2[, Sm])))
  Sm <- runValue(Rle(abf2[, Sm]))
  S0 <- list()
  i <- 1
  while (i < length(Sm)) {
    if(Sm[i + 1] < Sm[i] * (1 - minblockade)) {
      l0 <- which(Sm > Sm[i] * (1 - minblockade))
      l0 <- l0[l0 > i]
      if(length(l0) != 0) {
        Ei <- min(l0)
        S0 <- append(S0, list(c(end(gr[i]), start(gr[Ei]))))
        i <- Ei
      } else {
        i <- i + 1
      }
    } else {
      i <- i + 1
    }
  }
  
  Sigs <- lapply(S0, function(x) {
    L <- gr[end(gr) == x[1]]
    R <- gr[start(gr) == x[2]]
    gri <- gr[queryHits(findOverlaps(gr, IRanges(x[1], x[2]), type = "within"))]
    Main_ridge <- MainRidge(x = abf2[(x[1] + 1):(x[2] - 1), pA], bw = 0.005)
    SignalCurrent <- Main_ridge$MainPeak
    StageSD = abf2[(x[1] + 1):(x[2] - 1), sd(pA), Sm][, mean(V1)]
    SignalCurrentPercent <- Main_ridge$MainRidgeRatio * 100
    SignalCurrentWidth <- Main_ridge$MainRidgeWidth
    data.table(StartTime = abf2[x[1], Time], EndTime = abf2[x[2], Time], 
               LeftLength = width(L), RightLength = width(R), 
               BaseMean = mean2(c(abf2[start(L):end(L), pA], abf2[start(R):end(R), pA])),
               DeltaMean = abs(mean2(abf2[start(L):end(L), pA]) - mean2(abf2[start(R):end(R), pA])),
               StageSD = StageSD, 
               CurrentSD = abf2[(x[1] + 1):(x[2] - 1), sd(pA)], 
               Segments = sum(prop.table(width(gri)) > 0.1), 
               Valid = all(abf2[start(L):end(R), Valid]), 
               SignalCurrent = SignalCurrent, 
               SignalCurrentPercent = SignalCurrentPercent, 
               SignalCurrentWidth = SignalCurrentWidth)
  })
  Sigs <- do.call(rbind, Sigs)
  return(Sigs)
}

SignalCurrent <- function(x, abf) {
  target <- abf[, inrange(Time, x$StartTime, x$EndTime, incbounds = FALSE)]
  target <- IRanges(target)
  target <- as.data.table(target)
  target$ID <- seq_len(nrow(target))
  target <- target[, .(Loci = list(start:end)), ID]
  target <- target[, Loci]
  names(target) <- x$ID
  abf2 <- abf[do.call(c, target)]
  abf2$ID <- rep(names(target), mapply(length, target))
  abf2$BaseMean <- rep(x[, BaseMean], mapply(length, target))
  abf2[, Current := pA / BaseMean]
  abf2[, .(ID, Current)]
}


meta0 <- data.table(openxlsx::read.xlsx("./data/ChenShanchuan/20230612/数据记录20230612.xlsx", sheet = 2))
colnames(meta0) <- c("file_name", "start_time", "end_time", "amino_acid", "purpose", "priority", "sample", "base_line", "note", "temperature")
meta0 <- meta0[!is.na(amino_acid)]
meta0[, amino_acid := gsub(" ", "", amino_acid)]
meta0 <- meta0[!amino_acid %in% c("CPA", "AP")]
meta0 <- meta0[order(priority)]
meta0$abf_file <- mapply(meta0[, file_name], FUN = function(x) list.files("./data/ChenShanchuan", pattern = x, full.names = TRUE, recursive = TRUE))

meta <- meta0[, .(start_time = min(start_time), end_time = max(end_time)), abf_file]
meta[, file_name := gsub(".abf", "", basename(abf_file))]
setnames(meta, "abf_file", "file_path")


mclapply(1:nrow(meta), function(i) {
  print(i)
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"))) return()
  File <- meta[i, file_path]
  StartTime <- meta[i, start_time]
  EndTime <- meta[i, end_time]
  
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))) {
    abf <- readRDS(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  } else {
    abf <- tryCatch(readABF::readABF(File), error = function(e) NULL)
    if(is.null(abf)) return(NULL)
    abf <- as.data.table(as.data.frame(abf))
    if(ncol(abf) == 3) {
      colnames(abf) <- c("Time", "pA", "mV")
      abf[, Valid := round(mV) == 50 & pA > 0 & pA < 130]
    } else {
      colnames(abf) <- c("Time", "pA")
      abf[, Valid := pA > 0 & pA < 130]
    }
    abf <- abf[Time > StartTime * 60 & Time < EndTime * 60, ]
    abf$Sm <- KBins2(sig = abf$pA, pen.value = 1 - 1e-16, minseglen1 = 1, minseglen2 = 1)
    saveRDS(abf, paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  }
  
  # Raw Signal
  Sigs <- GetSignal(x = abf, minblockade = 0.1)
  Sigs[, DwellTime := EndTime - StartTime]
  Sigs <- Sigs[Valid == TRUE]
  if(nrow(Sigs) == 0) return(NULL)
  Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
  Sigs <- data.table(ID = paste(meta[i, file_name], Sigs[, seq_len(.N)], sep = "_"), Sigs)
  fwrite(Sigs, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"), sep = "\t", row.names = F, quote = F)
  
  # Signal Current
  Current <- SignalCurrent(Sigs, abf = abf)
  saveRDS(Current, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/SignalCurrent/SignalCurrent_", meta[i, file_name], ".Rds"))
  
  # Feature matrix
  Mat <- lapply(Sigs[, ID], function(j) {
    D <- density(Current[ID == j, Current], from = 0, to = 1, n = 1000, bw = 0.005)$y
    pA <- sort(Current[ID == j, Current])
    pA <- pA[round(seq(1, length(pA), length.out = 100))]
    c(D, pA)
  })
  Mat <- as.data.table(do.call(rbind, Mat))
  colnames(Mat) <- c(paste0("X", sprintf("%04d", seq_len(1000))), paste0("P", sprintf("%03d", 1:100)))
  Mat <- data.table(ID = Sigs[, ID], Mat)
  Mat <- merge(Mat, Sigs, by = "ID")
  setkey(Mat, ID)
  Mat <- Mat[Sigs$ID, ]
  fwrite(Mat, paste0("./analysis/21.ABFProcessing/03.MixtureAA/FeatureMatrix/FeatureMatrix_", meta[i, file_name], ".txt"), sep = "\t", quote = F)
}, mc.cores = 1)




meta0 <- data.table(openxlsx::read.xlsx("./data/ChenShanchuan/20230624/数据记录20230624.xlsx", sheet = 2))
colnames(meta0) <- c("file_name", "start_time", "end_time", "amino_acid", "purpose", "priority", "sample", "base_line", "note", "temperature")
meta0 <- meta0[!is.na(amino_acid)]
meta0[, amino_acid := gsub(" ", "", amino_acid)]
meta0 <- meta0[!amino_acid %in% c("CPA", "AP")]
meta0 <- meta0[order(priority)]
meta0$abf_file <- mapply(meta0[, file_name], FUN = function(x) list.files("./data/ChenShanchuan", pattern = x, full.names = TRUE, recursive = TRUE))

meta <- meta0[, .(start_time = min(start_time), end_time = max(end_time)), abf_file]
meta[, file_name := gsub(".abf", "", basename(abf_file))]
setnames(meta, "abf_file", "file_path")

meta <- meta[!file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", file_name, ".Rds"))]

mclapply(1:nrow(meta), function(i) {
  print(i)
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"))) return()
  File <- meta[i, file_path]
  StartTime <- meta[i, start_time]
  EndTime <- meta[i, end_time]
  
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))) {
    abf <- readRDS(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  } else {
    abf <- tryCatch(readABF::readABF(File), error = function(e) NULL)
    if(is.null(abf)) return(NULL)
    abf <- as.data.table(as.data.frame(abf))
    if(ncol(abf) == 3) {
      colnames(abf) <- c("Time", "pA", "mV")
      abf[, Valid := round(mV) == 50 & pA > 0 & pA < 130]
    } else {
      colnames(abf) <- c("Time", "pA")
      abf[, Valid := pA > 0 & pA < 130]
    }
    abf <- abf[Time > StartTime * 60 & Time < EndTime * 60, ]
    abf$Sm <- KBins2(sig = abf$pA, pen.value = 1 - 1e-16, minseglen1 = 1, minseglen2 = 1)
    saveRDS(abf, paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  }
  
  # Raw Signal
  Sigs <- GetSignal(x = abf, minblockade = 0.1)
  Sigs[, DwellTime := EndTime - StartTime]
  Sigs <- Sigs[Valid == TRUE]
  if(nrow(Sigs) == 0) return(NULL)
  Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
  Sigs <- data.table(ID = paste(meta[i, file_name], Sigs[, seq_len(.N)], sep = "_"), Sigs)
  fwrite(Sigs, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"), sep = "\t", row.names = F, quote = F)
  
  # Signal Current
  Current <- SignalCurrent(Sigs, abf = abf)
  saveRDS(Current, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/SignalCurrent/SignalCurrent_", meta[i, file_name], ".Rds"))
  
  # Feature matrix
  Mat <- lapply(Sigs[, ID], function(j) {
    D <- density(Current[ID == j, Current], from = 0, to = 1, n = 1000, bw = 0.005)$y
    pA <- sort(Current[ID == j, Current])
    pA <- pA[round(seq(1, length(pA), length.out = 100))]
    c(D, pA)
  })
  Mat <- as.data.table(do.call(rbind, Mat))
  colnames(Mat) <- c(paste0("X", sprintf("%04d", seq_len(1000))), paste0("P", sprintf("%03d", 1:100)))
  Mat <- data.table(ID = Sigs[, ID], Mat)
  Mat <- merge(Mat, Sigs, by = "ID")
  setkey(Mat, ID)
  Mat <- Mat[Sigs$ID, ]
  fwrite(Mat, paste0("./analysis/21.ABFProcessing/03.MixtureAA/FeatureMatrix/FeatureMatrix_", meta[i, file_name], ".txt"), sep = "\t", quote = F)
}, mc.cores = 1)






meta0 <- data.table(openxlsx::read.xlsx("./data/ChenShanchuan/20230701/数据记录20230701.xlsx", sheet = 2))
colnames(meta0) <- c("file_name", "start_time", "end_time", "amino_acid", "purpose", "priority", "sample", "base_line", "note", "temperature")
meta0 <- meta0[!is.na(amino_acid)]
meta0[, amino_acid := gsub(" ", "", amino_acid)]
meta0 <- meta0[!amino_acid %in% c("CPA", "AP")]
meta0 <- meta0[order(priority)]
meta0$abf_file <- mapply(meta0[, file_name], FUN = function(x) list.files("./data/ChenShanchuan", pattern = x, full.names = TRUE, recursive = TRUE))

meta <- meta0[, .(start_time = min(start_time), end_time = max(end_time)), abf_file]
meta[, file_name := gsub(".abf", "", basename(abf_file))]
setnames(meta, "abf_file", "file_path")

meta <- meta[!file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", file_name, ".Rds"))]

mclapply(1:nrow(meta), function(i) {
  print(i)
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"))) return()
  File <- meta[i, file_path]
  StartTime <- meta[i, start_time]
  EndTime <- meta[i, end_time]
  
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))) {
    abf <- readRDS(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  } else {
    abf <- tryCatch(readABF::readABF(File), error = function(e) NULL)
    if(is.null(abf)) return(NULL)
    abf <- as.data.table(as.data.frame(abf))
    if(ncol(abf) == 3) {
      colnames(abf) <- c("Time", "pA", "mV")
      abf[, Valid := round(mV) == 50 & pA > 0 & pA < 130]
    } else {
      colnames(abf) <- c("Time", "pA")
      abf[, Valid := pA > 0 & pA < 130]
    }
    abf <- abf[Time > StartTime * 60 & Time < EndTime * 60, ]
    abf$Sm <- KBins2(sig = abf$pA, pen.value = 1 - 1e-16, minseglen1 = 1, minseglen2 = 1)
    saveRDS(abf, paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  }
  
  # Raw Signal
  Sigs <- GetSignal(x = abf, minblockade = 0.1)
  Sigs[, DwellTime := EndTime - StartTime]
  Sigs <- Sigs[Valid == TRUE]
  if(nrow(Sigs) == 0) return(NULL)
  Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
  Sigs <- data.table(ID = paste(meta[i, file_name], Sigs[, seq_len(.N)], sep = "_"), Sigs)
  fwrite(Sigs, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"), sep = "\t", row.names = F, quote = F)
  
  # Signal Current
  Current <- SignalCurrent(Sigs, abf = abf)
  saveRDS(Current, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/SignalCurrent/SignalCurrent_", meta[i, file_name], ".Rds"))
  
  # Feature matrix
  Mat <- lapply(Sigs[, ID], function(j) {
    D <- density(Current[ID == j, Current], from = 0, to = 1, n = 1000, bw = 0.005)$y
    pA <- sort(Current[ID == j, Current])
    pA <- pA[round(seq(1, length(pA), length.out = 100))]
    c(D, pA)
  })
  Mat <- as.data.table(do.call(rbind, Mat))
  colnames(Mat) <- c(paste0("X", sprintf("%04d", seq_len(1000))), paste0("P", sprintf("%03d", 1:100)))
  Mat <- data.table(ID = Sigs[, ID], Mat)
  Mat <- merge(Mat, Sigs, by = "ID")
  setkey(Mat, ID)
  Mat <- Mat[Sigs$ID, ]
  fwrite(Mat, paste0("./analysis/21.ABFProcessing/03.MixtureAA/FeatureMatrix/FeatureMatrix_", meta[i, file_name], ".txt"), sep = "\t", quote = F)
}, mc.cores = 1)








meta0 <- data.table(openxlsx::read.xlsx("./data/ChenShanchuan/20230713/数据记录20230713.xlsx", sheet = 2))
colnames(meta0) <- c("file_name", "start_time", "end_time", "amino_acid", "purpose", "priority", "sample", "base_line", "note", "temperature")
meta0 <- meta0[!is.na(amino_acid)]
meta0[, amino_acid := gsub(" ", "", amino_acid)]
meta0 <- meta0[!amino_acid %in% c("CPA", "AP")]
meta0 <- meta0[order(priority)]
meta0$abf_file <- mapply(meta0[, file_name], FUN = function(x) list.files("./data/ChenShanchuan", pattern = x, full.names = TRUE, recursive = TRUE))

meta <- meta0[, .(start_time = min(start_time), end_time = max(end_time)), abf_file]
meta[, file_name := gsub(".abf", "", basename(abf_file))]
setnames(meta, "abf_file", "file_path")

meta <- meta[!file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", file_name, ".Rds"))]

mclapply(1:nrow(meta), function(i) {
  print(i)
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"))) return()
  File <- meta[i, file_path]
  StartTime <- meta[i, start_time]
  EndTime <- meta[i, end_time]
  
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))) {
    abf <- readRDS(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  } else {
    abf <- tryCatch(readABF::readABF(File), error = function(e) NULL)
    if(is.null(abf)) return(NULL)
    abf <- as.data.table(as.data.frame(abf))
    if(ncol(abf) == 3) {
      colnames(abf) <- c("Time", "pA", "mV")
      abf[, Valid := round(mV) == 50 & pA > 0 & pA < 130]
    } else {
      colnames(abf) <- c("Time", "pA")
      abf[, Valid := pA > 0 & pA < 130]
    }
    abf <- abf[Time > StartTime * 60 & Time < EndTime * 60, ]
    abf$Sm <- KBins2(sig = abf$pA, pen.value = 1 - 1e-16, minseglen1 = 1, minseglen2 = 1)
    saveRDS(abf, paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  }
  
  # Raw Signal
  Sigs <- GetSignal(x = abf, minblockade = 0.1)
  Sigs[, DwellTime := EndTime - StartTime]
  Sigs <- Sigs[Valid == TRUE]
  if(nrow(Sigs) == 0) return(NULL)
  Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
  Sigs <- data.table(ID = paste(meta[i, file_name], Sigs[, seq_len(.N)], sep = "_"), Sigs)
  fwrite(Sigs, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"), sep = "\t", row.names = F, quote = F)
  
  # Signal Current
  Current <- SignalCurrent(Sigs, abf = abf)
  saveRDS(Current, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/SignalCurrent/SignalCurrent_", meta[i, file_name], ".Rds"))
  
  # Feature matrix
  Mat <- lapply(Sigs[, ID], function(j) {
    D <- density(Current[ID == j, Current], from = 0, to = 1, n = 1000, bw = 0.005)$y
    pA <- sort(Current[ID == j, Current])
    pA <- pA[round(seq(1, length(pA), length.out = 100))]
    c(D, pA)
  })
  Mat <- as.data.table(do.call(rbind, Mat))
  colnames(Mat) <- c(paste0("X", sprintf("%04d", seq_len(1000))), paste0("P", sprintf("%03d", 1:100)))
  Mat <- data.table(ID = Sigs[, ID], Mat)
  Mat <- merge(Mat, Sigs, by = "ID")
  setkey(Mat, ID)
  Mat <- Mat[Sigs$ID, ]
  fwrite(Mat, paste0("./analysis/21.ABFProcessing/03.MixtureAA/FeatureMatrix/FeatureMatrix_", meta[i, file_name], ".txt"), sep = "\t", quote = F)
}, mc.cores = 1)








meta0 <- data.table(openxlsx::read.xlsx("./data/ChenShanchuan/20230715/数据记录20230715.xlsx", sheet = 2))
colnames(meta0) <- c("file_name", "start_time", "end_time", "amino_acid", "purpose", "priority", "sample", "base_line", "note", "temperature")
meta0 <- meta0[!is.na(amino_acid)]
meta0[, amino_acid := gsub(" ", "", amino_acid)]
meta0 <- meta0[!amino_acid %in% c("CPA", "AP")]
meta0 <- meta0[order(priority)]
meta0$abf_file <- mapply(meta0[, file_name], FUN = function(x) list.files("./data/ChenShanchuan", pattern = x, full.names = TRUE, recursive = TRUE))

meta <- meta0[, .(start_time = min(start_time), end_time = max(end_time)), abf_file]
meta[, file_name := gsub(".abf", "", basename(abf_file))]
setnames(meta, "abf_file", "file_path")

meta <- meta[!file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", file_name, ".Rds"))]

mclapply(1:nrow(meta), function(i) {
  print(i)
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"))) return()
  File <- meta[i, file_path]
  StartTime <- meta[i, start_time]
  EndTime <- meta[i, end_time]
  
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))) {
    abf <- readRDS(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  } else {
    abf <- tryCatch(readABF::readABF(File), error = function(e) NULL)
    if(is.null(abf)) return(NULL)
    abf <- as.data.table(as.data.frame(abf))
    if(ncol(abf) == 3) {
      colnames(abf) <- c("Time", "pA", "mV")
      abf[, Valid := round(mV) == 50 & pA > 0 & pA < 130]
    } else {
      colnames(abf) <- c("Time", "pA")
      abf[, Valid := pA > 0 & pA < 130]
    }
    abf <- abf[Time > StartTime * 60 & Time < EndTime * 60, ]
    abf$Sm <- KBins2(sig = abf$pA, pen.value = 1 - 1e-16, minseglen1 = 1, minseglen2 = 1)
    saveRDS(abf, paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  }
  
  # Raw Signal
  Sigs <- GetSignal(x = abf, minblockade = 0.1)
  Sigs[, DwellTime := EndTime - StartTime]
  Sigs <- Sigs[Valid == TRUE]
  if(nrow(Sigs) == 0) return(NULL)
  Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
  Sigs <- data.table(ID = paste(meta[i, file_name], Sigs[, seq_len(.N)], sep = "_"), Sigs)
  fwrite(Sigs, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"), sep = "\t", row.names = F, quote = F)
  
  # Signal Current
  Current <- SignalCurrent(Sigs, abf = abf)
  saveRDS(Current, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/SignalCurrent/SignalCurrent_", meta[i, file_name], ".Rds"))
  
  # Feature matrix
  Mat <- lapply(Sigs[, ID], function(j) {
    D <- density(Current[ID == j, Current], from = 0, to = 1, n = 1000, bw = 0.005)$y
    pA <- sort(Current[ID == j, Current])
    pA <- pA[round(seq(1, length(pA), length.out = 100))]
    c(D, pA)
  })
  Mat <- as.data.table(do.call(rbind, Mat))
  colnames(Mat) <- c(paste0("X", sprintf("%04d", seq_len(1000))), paste0("P", sprintf("%03d", 1:100)))
  Mat <- data.table(ID = Sigs[, ID], Mat)
  Mat <- merge(Mat, Sigs, by = "ID")
  setkey(Mat, ID)
  Mat <- Mat[Sigs$ID, ]
  fwrite(Mat, paste0("./analysis/21.ABFProcessing/03.MixtureAA/FeatureMatrix/FeatureMatrix_", meta[i, file_name], ".txt"), sep = "\t", quote = F)
}, mc.cores = 1)







meta0 <- data.table(openxlsx::read.xlsx("./data/ChenShanchuan/12种混合氨基酸检测/数据记录20230721.xlsx", sheet = 2))
colnames(meta0) <- c("file_name", "start_time", "end_time", "amino_acid", "purpose", "priority", "sample", "base_line", "note", "temperature")
meta0 <- meta0[file_name %in% gsub(".abf", "", list.files("data/ChenShanchuan/12种混合氨基酸检测", "abf"))]
meta0$abf_file <- mapply(meta0[, file_name], FUN = function(x) list.files("./data/ChenShanchuan", pattern = x, full.names = TRUE, recursive = TRUE))
meta <- meta0[, .(start_time = min(start_time), end_time = max(end_time)), abf_file]
meta[, file_name := gsub(".abf", "", basename(abf_file))]
setnames(meta, "abf_file", "file_path")
meta <- meta[!file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", file_name, ".Rds"))]

mclapply(1:nrow(meta), function(i) {
  print(i)
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"))) return()
  File <- meta[i, file_path]
  StartTime <- meta[i, start_time]
  EndTime <- meta[i, end_time]
  
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))) {
    abf <- readRDS(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  } else {
    abf <- tryCatch(readABF::readABF(File), error = function(e) NULL)
    if(is.null(abf)) return(NULL)
    abf <- as.data.table(as.data.frame(abf))
    if(ncol(abf) == 3) {
      colnames(abf) <- c("Time", "pA", "mV")
      abf[, Valid := round(mV) == 50 & pA > 0 & pA < 130]
    } else {
      colnames(abf) <- c("Time", "pA")
      abf[, Valid := pA > 0 & pA < 130]
    }
    abf <- abf[Time > StartTime * 60 & Time < EndTime * 60, ]
    abf$Sm <- KBins2(sig = abf$pA, pen.value = 1 - 1e-16, minseglen1 = 1, minseglen2 = 1)
    saveRDS(abf, paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  }
  
  # Raw Signal
  Sigs <- GetSignal(x = abf, minblockade = 0.1)
  Sigs[, DwellTime := EndTime - StartTime]
  Sigs <- Sigs[Valid == TRUE]
  if(nrow(Sigs) == 0) return(NULL)
  Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
  Sigs <- data.table(ID = paste(meta[i, file_name], Sigs[, seq_len(.N)], sep = "_"), Sigs)
  fwrite(Sigs, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"), sep = "\t", row.names = F, quote = F)
  
  # Signal Current
  Current <- SignalCurrent(Sigs, abf = abf)
  saveRDS(Current, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/SignalCurrent/SignalCurrent_", meta[i, file_name], ".Rds"))
  
  # Feature matrix
  Mat <- lapply(Sigs[, ID], function(j) {
    D <- density(Current[ID == j, Current], from = 0, to = 1, n = 1000, bw = 0.005)$y
    pA <- sort(Current[ID == j, Current])
    pA <- pA[round(seq(1, length(pA), length.out = 100))]
    c(D, pA)
  })
  Mat <- as.data.table(do.call(rbind, Mat))
  colnames(Mat) <- c(paste0("X", sprintf("%04d", seq_len(1000))), paste0("P", sprintf("%03d", 1:100)))
  Mat <- data.table(ID = Sigs[, ID], Mat)
  Mat <- merge(Mat, Sigs, by = "ID")
  setkey(Mat, ID)
  Mat <- Mat[Sigs$ID, ]
  fwrite(Mat, paste0("./analysis/21.ABFProcessing/03.MixtureAA/FeatureMatrix/FeatureMatrix_", meta[i, file_name], ".txt"), sep = "\t", quote = F)
}, mc.cores = 1)














meta0 <- data.table(openxlsx::read.xlsx("./data/ChenShanchuan/20230804/数据记录（0803）.xlsx", sheet = 1))
colnames(meta0) <- c("file_name", "start_time", "end_time", "amino_acid", "purpose", "sample", "base_line", "note")
meta0$abf_file <- mapply(meta0[, file_name], FUN = function(x) list.files("./data/ChenShanchuan", pattern = x, full.names = TRUE, recursive = TRUE))
meta <- meta0[, .(start_time = min(start_time), end_time = max(end_time)), abf_file]
meta[, file_name := gsub(".abf", "", basename(abf_file))]
setnames(meta, "abf_file", "file_path")
meta <- meta[!file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", file_name, ".Rds"))]

mclapply(1:nrow(meta), function(i) {
  print(i)
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"))) return()
  File <- meta[i, file_path]
  StartTime <- meta[i, start_time]
  EndTime <- meta[i, end_time]
  
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))) {
    abf <- readRDS(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  } else {
    abf <- tryCatch(readABF::readABF(File), error = function(e) NULL)
    if(is.null(abf)) return(NULL)
    abf <- as.data.table(as.data.frame(abf))
    if(ncol(abf) == 3) {
      colnames(abf) <- c("Time", "pA", "mV")
      abf[, Valid := round(mV) == 50 & pA > 0 & pA < 130]
    } else {
      colnames(abf) <- c("Time", "pA")
      abf[, Valid := pA > 0 & pA < 130]
    }
    abf <- abf[Time > StartTime * 60 & Time < EndTime * 60, ]
    abf$Sm <- KBins2(sig = abf$pA, pen.value = 1 - 1e-16, minseglen1 = 1, minseglen2 = 1)
    saveRDS(abf, paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  }
  
  # Raw Signal
  Sigs <- GetSignal(x = abf, minblockade = 0.1)
  Sigs[, DwellTime := EndTime - StartTime]
  Sigs <- Sigs[Valid == TRUE]
  if(nrow(Sigs) == 0) return(NULL)
  Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
  Sigs <- data.table(ID = paste(meta[i, file_name], Sigs[, seq_len(.N)], sep = "_"), Sigs)
  fwrite(Sigs, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"), sep = "\t", row.names = F, quote = F)
  
  # Signal Current
  Current <- SignalCurrent(Sigs, abf = abf)
  saveRDS(Current, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/SignalCurrent/SignalCurrent_", meta[i, file_name], ".Rds"))
  
  # Feature matrix
  Mat <- lapply(Sigs[, ID], function(j) {
    D <- density(Current[ID == j, Current], from = 0, to = 1, n = 1000, bw = 0.005)$y
    pA <- sort(Current[ID == j, Current])
    pA <- pA[round(seq(1, length(pA), length.out = 100))]
    c(D, pA)
  })
  Mat <- as.data.table(do.call(rbind, Mat))
  colnames(Mat) <- c(paste0("X", sprintf("%04d", seq_len(1000))), paste0("P", sprintf("%03d", 1:100)))
  Mat <- data.table(ID = Sigs[, ID], Mat)
  Mat <- merge(Mat, Sigs, by = "ID")
  setkey(Mat, ID)
  Mat <- Mat[Sigs$ID, ]
  fwrite(Mat, paste0("./analysis/21.ABFProcessing/03.MixtureAA/FeatureMatrix/FeatureMatrix_", meta[i, file_name], ".txt"), sep = "\t", quote = F)
}, mc.cores = 1)

















meta0 <- data.table(openxlsx::read.xlsx("./data/ChenShanchuan/20230816/副本数据记录20230816.xlsx", sheet = 1))
colnames(meta0) <- c("file_name", "start_time", "end_time", "amino_acid", "purpose", "level", "sample", "base_line", "note", "temp")
meta0 <- meta0[file_name %in% gsub(".abf", "", list.files("./data/ChenShanchuan/20230816", "abf"))]
meta0$abf_file <- mapply(meta0[, file_name], FUN = function(x) list.files("./data/ChenShanchuan/20230816", pattern = x, full.names = TRUE, recursive = TRUE))
meta <- meta0[, .(start_time = min(start_time), end_time = max(end_time)), abf_file]
meta[, file_name := gsub(".abf", "", basename(abf_file))]
setnames(meta, "abf_file", "file_path")
meta <- meta[!file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", file_name, ".Rds"))]

mclapply(1:nrow(meta), function(i) {
  print(i)
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"))) return()
  File <- meta[i, file_path]
  StartTime <- meta[i, start_time]
  EndTime <- meta[i, end_time]
  
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))) {
    abf <- readRDS(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  } else {
    abf <- tryCatch(readABF::readABF(File), error = function(e) NULL)
    if(is.null(abf)) return(NULL)
    abf <- as.data.table(as.data.frame(abf))
    if(ncol(abf) == 3) {
      colnames(abf) <- c("Time", "pA", "mV")
      abf[, Valid := round(mV) == 50 & pA > 0 & pA < 130]
    } else {
      colnames(abf) <- c("Time", "pA")
      abf[, Valid := pA > 0 & pA < 130]
    }
    abf <- abf[Time > StartTime * 60 & Time < EndTime * 60, ]
    abf$Sm <- KBins2(sig = abf$pA, pen.value = 1 - 1e-16, minseglen1 = 1, minseglen2 = 1)
    saveRDS(abf, paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  }
  
  # Raw Signal
  Sigs <- GetSignal(x = abf, minblockade = 0.1)
  Sigs[, DwellTime := EndTime - StartTime]
  Sigs <- Sigs[Valid == TRUE]
  if(nrow(Sigs) == 0) return(NULL)
  Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
  Sigs <- data.table(ID = paste(meta[i, file_name], Sigs[, seq_len(.N)], sep = "_"), Sigs)
  fwrite(Sigs, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"), sep = "\t", row.names = F, quote = F)
  
  # Signal Current
  Current <- SignalCurrent(Sigs, abf = abf)
  saveRDS(Current, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/SignalCurrent/SignalCurrent_", meta[i, file_name], ".Rds"))
  
  # Feature matrix
  Mat <- lapply(Sigs[, ID], function(j) {
    D <- density(Current[ID == j, Current], from = 0, to = 1, n = 1000, bw = 0.005)$y
    pA <- sort(Current[ID == j, Current])
    pA <- pA[round(seq(1, length(pA), length.out = 100))]
    c(D, pA)
  })
  Mat <- as.data.table(do.call(rbind, Mat))
  colnames(Mat) <- c(paste0("X", sprintf("%04d", seq_len(1000))), paste0("P", sprintf("%03d", 1:100)))
  Mat <- data.table(ID = Sigs[, ID], Mat)
  Mat <- merge(Mat, Sigs, by = "ID")
  setkey(Mat, ID)
  Mat <- Mat[Sigs$ID, ]
  fwrite(Mat, paste0("./analysis/21.ABFProcessing/03.MixtureAA/FeatureMatrix/FeatureMatrix_", meta[i, file_name], ".txt"), sep = "\t", quote = F)
}, mc.cores = 1)







meta0 <- data.table(openxlsx::read.xlsx("./data/ChenShanchuan/20230817/副本数据记录20230817.xlsx", sheet = 2))
colnames(meta0) <- c("file_name", "start_time", "end_time", "amino_acid", "purpose", "level", "sample", "base_line", "note", "temp")
meta0 <- meta0[file_name %in% gsub(".abf", "", list.files("./data/ChenShanchuan/20230817", "abf"))]
meta0$abf_file <- mapply(meta0[, file_name], FUN = function(x) list.files("./data/ChenShanchuan/20230817", pattern = x, full.names = TRUE, recursive = TRUE))
meta <- meta0[, .(start_time = min(start_time), end_time = max(end_time)), abf_file]
meta[, file_name := gsub(".abf", "", basename(abf_file))]
setnames(meta, "abf_file", "file_path")
meta <- meta[!file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", file_name, ".Rds"))]

mclapply(1:nrow(meta), function(i) {
  print(i)
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"))) return()
  File <- meta[i, file_path]
  StartTime <- meta[i, start_time]
  EndTime <- meta[i, end_time]
  
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))) {
    abf <- readRDS(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  } else {
    abf <- tryCatch(readABF::readABF(File), error = function(e) NULL)
    if(is.null(abf)) return(NULL)
    abf <- as.data.table(as.data.frame(abf))
    if(ncol(abf) == 3) {
      colnames(abf) <- c("Time", "pA", "mV")
      abf[, Valid := round(mV) == 50 & pA > 0 & pA < 130]
    } else {
      colnames(abf) <- c("Time", "pA")
      abf[, Valid := pA > 0 & pA < 130]
    }
    abf <- abf[Time > StartTime * 60 & Time < EndTime * 60, ]
    abf$Sm <- KBins2(sig = abf$pA, pen.value = 1 - 1e-16, minseglen1 = 1, minseglen2 = 1)
    saveRDS(abf, paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  }
  
  # Raw Signal
  Sigs <- GetSignal(x = abf, minblockade = 0.1)
  Sigs[, DwellTime := EndTime - StartTime]
  Sigs <- Sigs[Valid == TRUE]
  if(nrow(Sigs) == 0) return(NULL)
  Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
  Sigs <- data.table(ID = paste(meta[i, file_name], Sigs[, seq_len(.N)], sep = "_"), Sigs)
  fwrite(Sigs, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"), sep = "\t", row.names = F, quote = F)
  
  # Signal Current
  Current <- SignalCurrent(Sigs, abf = abf)
  saveRDS(Current, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/SignalCurrent/SignalCurrent_", meta[i, file_name], ".Rds"))
  
  # Feature matrix
  Mat <- lapply(Sigs[, ID], function(j) {
    D <- density(Current[ID == j, Current], from = 0, to = 1, n = 1000, bw = 0.005)$y
    pA <- sort(Current[ID == j, Current])
    pA <- pA[round(seq(1, length(pA), length.out = 100))]
    c(D, pA)
  })
  Mat <- as.data.table(do.call(rbind, Mat))
  colnames(Mat) <- c(paste0("X", sprintf("%04d", seq_len(1000))), paste0("P", sprintf("%03d", 1:100)))
  Mat <- data.table(ID = Sigs[, ID], Mat)
  Mat <- merge(Mat, Sigs, by = "ID")
  setkey(Mat, ID)
  Mat <- Mat[Sigs$ID, ]
  fwrite(Mat, paste0("./analysis/21.ABFProcessing/03.MixtureAA/FeatureMatrix/FeatureMatrix_", meta[i, file_name], ".txt"), sep = "\t", quote = F)
}, mc.cores = 1)








meta0 <- data.table(openxlsx::read.xlsx("./data/ChenShanchuan/20230817_2/副本数据记录20230817新.xlsx", sheet = 2))
colnames(meta0) <- c("file_name", "start_time", "end_time", "amino_acid", "purpose", "level", "sample", "base_line", "note", "temp")
meta0 <- meta0[file_name %in% gsub(".abf", "", list.files("./data/ChenShanchuan/20230817_2", "abf"))]
meta0$abf_file <- mapply(meta0[, file_name], FUN = function(x) list.files("./data/ChenShanchuan/20230817_2", pattern = x, full.names = TRUE, recursive = TRUE))
meta <- meta0[, .(start_time = min(start_time), end_time = max(end_time)), abf_file]
meta[, file_name := gsub(".abf", "", basename(abf_file))]
setnames(meta, "abf_file", "file_path")
meta <- meta[!file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", file_name, ".Rds"))]

mclapply(1:nrow(meta), function(i) {
  print(i)
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"))) return()
  File <- meta[i, file_path]
  StartTime <- meta[i, start_time]
  EndTime <- meta[i, end_time]
  
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))) {
    abf <- readRDS(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  } else {
    abf <- tryCatch(readABF::readABF(File), error = function(e) NULL)
    if(is.null(abf)) return(NULL)
    abf <- as.data.table(as.data.frame(abf))
    if(ncol(abf) == 3) {
      colnames(abf) <- c("Time", "pA", "mV")
      abf[, Valid := round(mV) == 50 & pA > 0 & pA < 130]
    } else {
      colnames(abf) <- c("Time", "pA")
      abf[, Valid := pA > 0 & pA < 130]
    }
    abf <- abf[Time > StartTime * 60 & Time < EndTime * 60, ]
    abf$Sm <- KBins2(sig = abf$pA, pen.value = 1 - 1e-16, minseglen1 = 1, minseglen2 = 1)
    saveRDS(abf, paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  }
  
  # Raw Signal
  Sigs <- GetSignal(x = abf, minblockade = 0.1)
  Sigs[, DwellTime := EndTime - StartTime]
  Sigs <- Sigs[Valid == TRUE]
  if(nrow(Sigs) == 0) return(NULL)
  Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
  Sigs <- data.table(ID = paste(meta[i, file_name], Sigs[, seq_len(.N)], sep = "_"), Sigs)
  fwrite(Sigs, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"), sep = "\t", row.names = F, quote = F)
  
  # Signal Current
  Current <- SignalCurrent(Sigs, abf = abf)
  saveRDS(Current, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/SignalCurrent/SignalCurrent_", meta[i, file_name], ".Rds"))
  
  # Feature matrix
  Mat <- lapply(Sigs[, ID], function(j) {
    D <- density(Current[ID == j, Current], from = 0, to = 1, n = 1000, bw = 0.005)$y
    pA <- sort(Current[ID == j, Current])
    pA <- pA[round(seq(1, length(pA), length.out = 100))]
    c(D, pA)
  })
  Mat <- as.data.table(do.call(rbind, Mat))
  colnames(Mat) <- c(paste0("X", sprintf("%04d", seq_len(1000))), paste0("P", sprintf("%03d", 1:100)))
  Mat <- data.table(ID = Sigs[, ID], Mat)
  Mat <- merge(Mat, Sigs, by = "ID")
  setkey(Mat, ID)
  Mat <- Mat[Sigs$ID, ]
  fwrite(Mat, paste0("./analysis/21.ABFProcessing/03.MixtureAA/FeatureMatrix/FeatureMatrix_", meta[i, file_name], ".txt"), sep = "\t", quote = F)
}, mc.cores = 1)















meta0 <- data.table(openxlsx::read.xlsx("./data/ChenShanchuan/20230818/副本数据记录20230818.xlsx", sheet = 2))
colnames(meta0) <- c("file_name", "start_time", "end_time", "amino_acid", "purpose", "level", "sample", "base_line", "note", "temp")
meta0 <- meta0[file_name %in% gsub(".abf", "", list.files("./data/ChenShanchuan/20230818", "abf"))]
meta0$abf_file <- mapply(meta0[, file_name], FUN = function(x) list.files("./data/ChenShanchuan/20230818", pattern = x, full.names = TRUE, recursive = TRUE))
meta <- meta0[, .(start_time = min(start_time), end_time = max(end_time)), abf_file]
meta[, file_name := gsub(".abf", "", basename(abf_file))]
setnames(meta, "abf_file", "file_path")
meta <- meta[!file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", file_name, ".Rds"))]

mclapply(1:nrow(meta), function(i) {
  print(i)
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"))) return()
  File <- meta[i, file_path]
  StartTime <- meta[i, start_time]
  EndTime <- meta[i, end_time]
  
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))) {
    abf <- readRDS(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  } else {
    abf <- tryCatch(readABF::readABF(File), error = function(e) NULL)
    if(is.null(abf)) return(NULL)
    abf <- as.data.table(as.data.frame(abf))
    if(ncol(abf) == 3) {
      colnames(abf) <- c("Time", "pA", "mV")
      abf[, Valid := round(mV) == 50 & pA > 0 & pA < 130]
    } else {
      colnames(abf) <- c("Time", "pA")
      abf[, Valid := pA > 0 & pA < 130]
    }
    abf <- abf[Time > StartTime * 60 & Time < EndTime * 60, ]
    abf$Sm <- KBins2(sig = abf$pA, pen.value = 1 - 1e-16, minseglen1 = 1, minseglen2 = 1)
    saveRDS(abf, paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  }
  
  # Raw Signal
  Sigs <- GetSignal(x = abf, minblockade = 0.1)
  Sigs[, DwellTime := EndTime - StartTime]
  Sigs <- Sigs[Valid == TRUE]
  if(nrow(Sigs) == 0) return(NULL)
  Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
  Sigs <- data.table(ID = paste(meta[i, file_name], Sigs[, seq_len(.N)], sep = "_"), Sigs)
  fwrite(Sigs, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"), sep = "\t", row.names = F, quote = F)
  
  # Signal Current
  Current <- SignalCurrent(Sigs, abf = abf)
  saveRDS(Current, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/SignalCurrent/SignalCurrent_", meta[i, file_name], ".Rds"))
  
  # Feature matrix
  Mat <- lapply(Sigs[, ID], function(j) {
    D <- density(Current[ID == j, Current], from = 0, to = 1, n = 1000, bw = 0.005)$y
    pA <- sort(Current[ID == j, Current])
    pA <- pA[round(seq(1, length(pA), length.out = 100))]
    c(D, pA)
  })
  Mat <- as.data.table(do.call(rbind, Mat))
  colnames(Mat) <- c(paste0("X", sprintf("%04d", seq_len(1000))), paste0("P", sprintf("%03d", 1:100)))
  Mat <- data.table(ID = Sigs[, ID], Mat)
  Mat <- merge(Mat, Sigs, by = "ID")
  setkey(Mat, ID)
  Mat <- Mat[Sigs$ID, ]
  fwrite(Mat, paste0("./analysis/21.ABFProcessing/03.MixtureAA/FeatureMatrix/FeatureMatrix_", meta[i, file_name], ".txt"), sep = "\t", quote = F)
}, mc.cores = 1)















meta0 <- data.table(openxlsx::read.xlsx("./data/ChenShanchuan/20230822/副本数据记录20230822.xlsx", sheet = 2))
colnames(meta0) <- c("file_name", "start_time", "end_time", "amino_acid", "purpose", "level", "sample", "base_line", "note", "temp")
meta0 <- meta0[file_name %in% gsub(".abf", "", list.files("./data/ChenShanchuan/20230822", "abf"))]
meta0$abf_file <- mapply(meta0[, file_name], FUN = function(x) list.files("./data/ChenShanchuan/20230822", pattern = x, full.names = TRUE, recursive = TRUE))
meta <- meta0[, .(start_time = min(start_time), end_time = max(end_time)), abf_file]
meta[, file_name := gsub(".abf", "", basename(abf_file))]
setnames(meta, "abf_file", "file_path")
meta <- meta[!file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", file_name, ".Rds"))]

mclapply(1:nrow(meta), function(i) {
  print(i)
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"))) return()
  File <- meta[i, file_path]
  StartTime <- meta[i, start_time]
  EndTime <- meta[i, end_time]
  
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))) {
    abf <- readRDS(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  } else {
    abf <- tryCatch(readABF::readABF(File), error = function(e) NULL)
    if(is.null(abf)) return(NULL)
    abf <- as.data.table(as.data.frame(abf))
    if(ncol(abf) == 3) {
      colnames(abf) <- c("Time", "pA", "mV")
      abf[, Valid := round(mV) == 50 & pA > 0 & pA < 130]
    } else {
      colnames(abf) <- c("Time", "pA")
      abf[, Valid := pA > 0 & pA < 130]
    }
    abf <- abf[Time > StartTime * 60 & Time < EndTime * 60, ]
    abf$Sm <- KBins2(sig = abf$pA, pen.value = 1 - 1e-16, minseglen1 = 1, minseglen2 = 1)
    saveRDS(abf, paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  }
  
  # Raw Signal
  Sigs <- GetSignal(x = abf, minblockade = 0.1)
  Sigs[, DwellTime := EndTime - StartTime]
  Sigs <- Sigs[Valid == TRUE]
  if(nrow(Sigs) == 0) return(NULL)
  Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
  Sigs <- data.table(ID = paste(meta[i, file_name], Sigs[, seq_len(.N)], sep = "_"), Sigs)
  fwrite(Sigs, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"), sep = "\t", row.names = F, quote = F)
  
  # Signal Current
  Current <- SignalCurrent(Sigs, abf = abf)
  saveRDS(Current, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/SignalCurrent/SignalCurrent_", meta[i, file_name], ".Rds"))
  
  # Feature matrix
  Mat <- lapply(Sigs[, ID], function(j) {
    D <- density(Current[ID == j, Current], from = 0, to = 1, n = 1000, bw = 0.005)$y
    pA <- sort(Current[ID == j, Current])
    pA <- pA[round(seq(1, length(pA), length.out = 100))]
    c(D, pA)
  })
  Mat <- as.data.table(do.call(rbind, Mat))
  colnames(Mat) <- c(paste0("X", sprintf("%04d", seq_len(1000))), paste0("P", sprintf("%03d", 1:100)))
  Mat <- data.table(ID = Sigs[, ID], Mat)
  Mat <- merge(Mat, Sigs, by = "ID")
  setkey(Mat, ID)
  Mat <- Mat[Sigs$ID, ]
  fwrite(Mat, paste0("./analysis/21.ABFProcessing/03.MixtureAA/FeatureMatrix/FeatureMatrix_", meta[i, file_name], ".txt"), sep = "\t", quote = F)
}, mc.cores = 1)


meta0 <- data.table(openxlsx::read.xlsx("./data/ChenShanchuan/20230822/副本数据记录20230822.xlsx", sheet = 2))
colnames(meta0) <- c("file_name", "start_time", "end_time", "amino_acid", "purpose", "level", "sample", "base_line", "note", "temp")
meta0 <- meta0[file_name %in% gsub(".abf", "", list.files("./data/ChenShanchuan/20230822", "abf"))]
meta0$abf_file <- mapply(meta0[, file_name], FUN = function(x) list.files("./data/ChenShanchuan/20230822", pattern = x, full.names = TRUE, recursive = TRUE))
meta <- meta0[, .(start_time = min(start_time), end_time = max(end_time)), abf_file]
meta[, file_name := gsub(".abf", "", basename(abf_file))]
setnames(meta, "abf_file", "file_path")
meta <- meta[4, ]
meta[, file_name := "20230822_0006_2"]
meta[, start_time := 12]


mclapply(1:nrow(meta), function(i) {
  print(i)
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"))) return()
  File <- meta[i, file_path]
  StartTime <- meta[i, start_time]
  EndTime <- meta[i, end_time]
  
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))) {
    abf <- readRDS(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  } else {
    abf <- tryCatch(readABF::readABF(File), error = function(e) NULL)
    if(is.null(abf)) return(NULL)
    abf <- as.data.table(as.data.frame(abf))
    if(ncol(abf) == 3) {
      colnames(abf) <- c("Time", "pA", "mV")
      abf[, Valid := round(mV) == 50 & pA > 0 & pA < 130]
    } else {
      colnames(abf) <- c("Time", "pA")
      abf[, Valid := pA > 0 & pA < 130]
    }
    abf <- abf[Time > StartTime * 60 & Time < EndTime * 60, ]
    abf$Sm <- KBins2(sig = abf$pA, pen.value = 1 - 1e-16, minseglen1 = 1, minseglen2 = 1)
    saveRDS(abf, paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  }
  
  # Raw Signal
  Sigs <- GetSignal(x = abf, minblockade = 0.1)
  Sigs[, DwellTime := EndTime - StartTime]
  Sigs <- Sigs[Valid == TRUE]
  if(nrow(Sigs) == 0) return(NULL)
  Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
  Sigs <- data.table(ID = paste(meta[i, file_name], Sigs[, seq_len(.N)], sep = "_"), Sigs)
  fwrite(Sigs, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"), sep = "\t", row.names = F, quote = F)
  
  # Signal Current
  Current <- SignalCurrent(Sigs, abf = abf)
  saveRDS(Current, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/SignalCurrent/SignalCurrent_", meta[i, file_name], ".Rds"))
  
  # Feature matrix
  Mat <- lapply(Sigs[, ID], function(j) {
    D <- density(Current[ID == j, Current], from = 0, to = 1, n = 1000, bw = 0.005)$y
    pA <- sort(Current[ID == j, Current])
    pA <- pA[round(seq(1, length(pA), length.out = 100))]
    c(D, pA)
  })
  Mat <- as.data.table(do.call(rbind, Mat))
  colnames(Mat) <- c(paste0("X", sprintf("%04d", seq_len(1000))), paste0("P", sprintf("%03d", 1:100)))
  Mat <- data.table(ID = Sigs[, ID], Mat)
  Mat <- merge(Mat, Sigs, by = "ID")
  setkey(Mat, ID)
  Mat <- Mat[Sigs$ID, ]
  fwrite(Mat, paste0("./analysis/21.ABFProcessing/03.MixtureAA/FeatureMatrix/FeatureMatrix_", meta[i, file_name], ".txt"), sep = "\t", quote = F)
}, mc.cores = 1)



















meta0 <- data.table(openxlsx::read.xlsx("./data/ChenShanchuan/20230823/副本数据记录20230823.xlsx", sheet = 2))
colnames(meta0) <- c("file_name", "start_time", "end_time", "amino_acid", "purpose", "level", "sample", "base_line", "note", "temp")
meta0 <- meta0[file_name %in% gsub(".abf", "", list.files("./data/ChenShanchuan/20230823", "abf"))]
meta0$abf_file <- mapply(meta0[, file_name], FUN = function(x) list.files("./data/ChenShanchuan/20230823", pattern = x, full.names = TRUE, recursive = TRUE))
meta <- meta0[, .(start_time = min(start_time), end_time = max(end_time)), abf_file]
meta[, file_name := gsub(".abf", "", basename(abf_file))]
setnames(meta, "abf_file", "file_path")
meta <- meta[!file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", file_name, ".Rds"))]

mclapply(1:nrow(meta), function(i) {
  print(i)
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"))) return()
  File <- meta[i, file_path]
  StartTime <- meta[i, start_time]
  EndTime <- meta[i, end_time]
  
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))) {
    abf <- readRDS(paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  } else {
    abf <- tryCatch(readABF::readABF(File), error = function(e) NULL)
    if(is.null(abf)) return(NULL)
    abf <- as.data.table(as.data.frame(abf))
    if(ncol(abf) == 3) {
      colnames(abf) <- c("Time", "pA", "mV")
      abf[, Valid := round(mV) == 50 & pA > 0 & pA < 130]
    } else {
      colnames(abf) <- c("Time", "pA")
      abf[, Valid := pA > 0 & pA < 130]
    }
    abf <- abf[Time > StartTime * 60 & Time < EndTime * 60, ]
    abf$Sm <- KBins2(sig = abf$pA, pen.value = 1 - 1e-16, minseglen1 = 1, minseglen2 = 1)
    saveRDS(abf, paste0("./analysis/21.ABFProcessing/03.MixtureAA/ABF/ABF_", meta[i, file_name], ".Rds"))
  }
  
  # Raw Signal
  Sigs <- GetSignal(x = abf, minblockade = 0.1)
  Sigs[, DwellTime := EndTime - StartTime]
  Sigs <- Sigs[Valid == TRUE]
  if(nrow(Sigs) == 0) return(NULL)
  Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
  Sigs <- data.table(ID = paste(meta[i, file_name], Sigs[, seq_len(.N)], sep = "_"), Sigs)
  fwrite(Sigs, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"), sep = "\t", row.names = F, quote = F)
  
  # Signal Current
  Current <- SignalCurrent(Sigs, abf = abf)
  saveRDS(Current, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/SignalCurrent/SignalCurrent_", meta[i, file_name], ".Rds"))
  
  # Feature matrix
  Mat <- lapply(Sigs[, ID], function(j) {
    D <- density(Current[ID == j, Current], from = 0, to = 1, n = 1000, bw = 0.005)$y
    pA <- sort(Current[ID == j, Current])
    pA <- pA[round(seq(1, length(pA), length.out = 100))]
    c(D, pA)
  })
  Mat <- as.data.table(do.call(rbind, Mat))
  colnames(Mat) <- c(paste0("X", sprintf("%04d", seq_len(1000))), paste0("P", sprintf("%03d", 1:100)))
  Mat <- data.table(ID = Sigs[, ID], Mat)
  Mat <- merge(Mat, Sigs, by = "ID")
  setkey(Mat, ID)
  Mat <- Mat[Sigs$ID, ]
  fwrite(Mat, paste0("./analysis/21.ABFProcessing/03.MixtureAA/FeatureMatrix/FeatureMatrix_", meta[i, file_name], ".txt"), sep = "\t", quote = F)
}, mc.cores = 1)

