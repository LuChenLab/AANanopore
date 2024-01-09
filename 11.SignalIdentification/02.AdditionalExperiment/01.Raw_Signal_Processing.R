setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")

library(data.table)
library(Biostrings)
library(IRanges)
library(ggplot2)
library(parallel)
library(changepoint)

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
    SignalCurrent <- with(density(abf2[(x[1] + 1):(x[2] - 1), pA], n = 10000, adjust = 1), x[which.max(y)])
    SignalCurrent.5 <- with(density(abf2[(x[1] + 1):(x[2] - 1), pA], n = 10000, adjust = .75), x[which.max(y)])
    L2Ratio = mean(abf2[(x[1] + 1):(x[2] - 1), Sm] < SignalCurrent * (1 - minblockade))
    StageSD = abf2[(x[1] + 1):(x[2] - 1), sd(pA), Sm][, mean(V1)]
    SignalCurrentPercent <- mean(abs(abf2[(x[1] + 1):(x[2] - 1), Sm] - SignalCurrent) < StageSD) * 100
    data.table(StartTime = abf2[x[1], Time], EndTime = abf2[x[2], Time], 
               LeftLength = width(L), RightLength = width(R), 
               BaseMean = mean2(c(abf2[start(L):end(L), pA], abf2[start(R):end(R), pA])),
               DeltaMean = abs(mean2(abf2[start(L):end(L), pA]) - mean2(abf2[start(R):end(R), pA])),
               StageSD = StageSD, 
               CurrentSD = abf2[(x[1] + 1):(x[2] - 1), sd(pA)], 
               Segments = sum(prop.table(width(gri)) > 0.1), Valid = all(abf2[start(L):end(R), Valid]), 
               SignalCurrent = SignalCurrent, SignalCurrent.5 = SignalCurrent.5, L2Ratio = L2Ratio, 
               SignalCurrentPercent = SignalCurrentPercent)
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


meta0 <- data.table(openxlsx::read.xlsx("./data/ChenShanchuan/20230521/副本副本数据记录20230521新.xlsx", sheet = 2))
colnames(meta0) <- c("file_name", "start_time", "end_time", "amino_acid", "purpose", "priority", "sample", "base_line", "note", "temperature")
meta0 <- meta0[!is.na(amino_acid)]
meta0[, amino_acid := gsub(" ", "", amino_acid)]
meta0 <- meta0[!amino_acid %in% c("CPA", "AP")]
meta0 <- meta0[order(priority)]
meta0$abf_file <- mapply(meta0[, file_name], FUN = function(x) list.files("./data/ChenShanchuan", pattern = x, full.names = TRUE, recursive = TRUE))

meta <- meta0[, .(start_time = min(start_time), end_time = max(end_time)), abf_file]
meta[, file_name := gsub(".abf", "", basename(abf_file))]

mclapply(1:nrow(meta), function(i) {
  print(i)
  if(file.exists(paste0("./analysis/11.SignalIdentification/FeatureMatrix/FeatureMatrix_", meta[i, file_name], ".txt"))) {
    return(NULL)
  }
  File <- meta[i, abf_file]
  StartTime <- meta[i, start_time]
  EndTime <- meta[i, end_time]
  
  if(file.exists(paste0("./analysis/11.SignalIdentification/ABF/ABF_", meta[i, file_name], ".Rds"))) {
    abf <- readRDS(paste0("./analysis/11.SignalIdentification/ABF/ABF_", meta[i, file_name], ".Rds"))
  } else {
    abf <- readABF::readABF(File)
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
  }
  saveRDS(abf, file = paste0("./analysis/11.SignalIdentification/ABF/ABF_", meta[i, file_name], ".Rds"))
  Sigs <- GetSignal(x = abf)
  Sigs[, DwellTime := EndTime - StartTime]
  Sigs <- Sigs[Valid == TRUE]
  if(nrow(Sigs) == 0) return(NULL)
  Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
  Sigs <- data.table(ID = paste(meta[i, file_name], Sigs[, seq_len(.N)], sep = "_"), Sigs)
  Current <- SignalCurrent(Sigs, abf = abf)
  saveRDS(Current, file = paste0("./analysis/11.SignalIdentification/SignalCurrent/SignalCurrent_", meta[i, file_name], ".Rds"))
  fwrite(Sigs, file = paste0("./analysis/11.SignalIdentification/Signal/RawSignal_", meta[i, file_name], ".txt"), sep = "\t", row.names = F, quote = F)
  
  # Feature matrix
  Mat <- lapply(Sigs[, ID], function(j) {
    D <- density(Current[ID == j, Current], from = 0, to = 1, n = 500, adjust = 0.5)$y
    pA <- sort(Current[ID == j, Current])
    pA <- pA[round(seq(1, length(pA), length.out = 100))]
    c(D, pA)
  })
  Mat <- as.data.table(do.call(rbind, Mat))
  colnames(Mat) <- c(paste0("X", sprintf("%03d", 1:500)), paste0("P", sprintf("%03d", 1:100)))
  Mat <- data.table(ID = Sigs[, ID], Mat)
  Mat <- merge(Mat, Sigs, by = "ID")
  setkey(Mat, ID)
  Mat <- Mat[Sigs$ID, ]
  fwrite(Mat, paste0("./analysis/11.SignalIdentification/FeatureMatrix/FeatureMatrix_", meta[i, file_name], ".txt"), sep = "\t", quote = F)
}, mc.cores = 1)



