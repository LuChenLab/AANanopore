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

GetSignal <- function(abf, minblockade = 0.1, cores = 1) {
  fun1 <- function(x, y, minblockade = 0.1, cores = cores) {
    y2 <- Rle(y)
    uR <- IRanges(diff(runValue(y2)) >= 0)
    dR <- IRanges(diff(runValue(y2)) <= 0)
    if(start(uR[1]) == 1) uR <- uR[-1]
    if(max(end(dR)) > max(end(uR))) dR <- dR[-c(length(dR))]
    stopifnot(length(uR) == length(dR))
    
    res <- IRanges(end = end(uR), start = start(dR))
    res <- IRanges(start = mcmapply(function(x) sum(runLength(y2)[1:x]), start(res), mc.cores = cores), 
                   end = mcmapply(function(x) sum(runLength(y2)[1:x]), end(res), mc.cores = cores), 
                   yL = runValue(y2)[start(res)], 
                   yR = runValue(y2)[end(res) + 1])
    ISminblockade <- mcols(res)$yR >= mcols(res)$yL * (1 - minblockade) & mcols(res)$yL >= mcols(res)$yR * (1 - minblockade)
    if(any(!ISminblockade)) {
      res1 <- res[ISminblockade]
      res2 <- res[!ISminblockade]
      mcols(res2)$end2 <- mcmapply(FUN = function(i) {
        ys <- start(IRanges(y >= mcols(res2)$yL[i] * (1 - minblockade)))
        min(ys[ys > start(res2[i])]) - 1
      }, seq_along(res2), mc.cores = cores)
      
      mcols(res2)$start2 <- mcmapply(FUN = function(i) {
        ys <- end(IRanges(y >= mcols(res2)$yR[i] * (1 - minblockade)))
        max(ys[ys < end(res2[i])])
      }, seq_along(res2), mc.cores = cores)
      
      res2 <- unique(c(IRanges(start = start(res2), end = ifelse(is.infinite(mcols(res2)$end2), length(x), mcols(res2)$end2)), 
                       IRanges(end = end(res2), start = ifelse(is.infinite(mcols(res2)$start2), 0, mcols(res2)$start2))))
      res2 <- res2[start(res2) != 0 & end(res2) != length(x)]
      res2 <- IRanges(res2, yL = y[start(res2)], yR = y[end(res2) + 1])
      res2 <- res2[mcols(res2)$yR >= mcols(res2)$yL * (1 - minblockade) & mcols(res2)$yL >= mcols(res2)$yR * (1 - minblockade)]
      res <- sort(c(res1, res2))
    }
    # data.table(start = x[start(res)], end = x[end(res)], as.data.frame(mcols(res)))
    return(res)
  }
  
  gr <- IRanges(end = cumsum(runLength(Rle(abf[, Sm]))), width = runLength(Rle(abf[, Sm])))
  S0 <- abf[, suppressWarnings(fun1(Time, Sm, minblockade = minblockade, cores = cores))]
  
  Sigs <- mclapply(seq_along(S0), function(j) {
    # print(j)
    x <- c(start(S0[j]), end(S0[j]) + 1)
    L <- gr[end(gr) == x[1]]
    R <- gr[start(gr) == x[2]]
    BaseMean <- mean2(c(abf[start(L):end(L), pA], abf[start(R):end(R), pA]))
    StartTime <- suppressWarnings(abf[x[1]:x[2]][Sm < BaseMean * (1 - minblockade), min(Time)])
    EndTime <- suppressWarnings(abf[x[1]:x[2]][Sm < BaseMean * (1 - minblockade), max(Time)])
    if(any(is.infinite(c(StartTime, EndTime)))) return(NULL)
    gri <- gr[queryHits(findOverlaps(gr, IRanges(x[1], x[2]), type = "within"))]
    Main_ridge <- MainRidge(x = abf[between(Time, StartTime, EndTime), pA], bw = 0.005)
    SignalCurrent <- Main_ridge$MainPeak
    StageSD = abf[between(Time, StartTime, EndTime), sd(pA), Sm][, mean(V1)]
    SignalCurrentPercent <- Main_ridge$MainRidgeRatio * 100
    SignalCurrentWidth <- Main_ridge$MainRidgeWidth
    data.table(StartTime = StartTime, EndTime = EndTime, 
               LeftLength = width(L), RightLength = width(R), 
               BaseMean = BaseMean,
               DeltaMean = abs(mean2(abf[start(L):end(L), pA]) - mean2(abf[start(R):end(R), pA])),
               StageSD = StageSD, 
               CurrentSD = abf[between(Time, StartTime, EndTime), sd(pA)], 
               Segments = max(1, sum(prop.table(width(gri)) > mean(prop.table(width(gri))))), 
               Valid = all(abf[start(L):end(R), Valid]), 
               SignalCurrent = SignalCurrent, 
               SignalCurrentPercent = SignalCurrentPercent, 
               SignalCurrentWidth = SignalCurrentWidth)
  }, mc.cores = cores)
  Sigs <- do.call(rbind, Sigs)
  Sigs[, DwellTime := EndTime - StartTime]
  Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
  return(Sigs)
}

SignalCurrent <- function(x, abf, cores = 1) {
  target <- abf[inrange(Time, x$StartTime, x$EndTime, incbounds = FALSE)]
  target <- mclapply(seq_len(nrow(x)), function(i) {
    j <- target[between(Time, x[i, StartTime], x[i, EndTime])]
    j[, .(ID = x[i, ID], Current = pA / x[i, BaseMean])]
  }, mc.cores = cores)
  do.call(rbind, target)
}


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


lapply(1:nrow(meta), function(i) {
  print(i)
  if(file.exists(paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"))) return(NULL)
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
  Sigs <- GetSignal(abf = abf, minblockade = 0.1, cores = 20)
  Sigs[, DwellTime := EndTime - StartTime]
  Sigs <- Sigs[Valid == TRUE]
  if(nrow(Sigs) == 0) return(NULL)
  Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
  Sigs <- data.table(ID = paste(meta[i, file_name], Sigs[, seq_len(.N)], sep = "_"), Sigs)
  fwrite(Sigs, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/RawSignal/RawSignal_", meta[i, file_name], ".txt"), sep = "\t", row.names = F, quote = F)
  
  # Signal Current
  Current <- SignalCurrent(Sigs, abf = abf, cores = 20)
  saveRDS(Current, file = paste0("./analysis/21.ABFProcessing/03.MixtureAA/SignalCurrent/SignalCurrent_", meta[i, file_name], ".Rds"))
  
  # Feature matrix
  Mat <- mclapply(Sigs[, ID], function(j) {
    D <- density(Current[ID == j, Current], from = 0, to = 1, n = 1000, bw = 0.005)$y
    pA <- sort(Current[ID == j, Current])
    pA <- pA[round(seq(1, length(pA), length.out = 100))]
    c(D, pA)
  }, mc.cores = 20)
  Mat <- as.data.table(do.call(rbind, Mat))
  colnames(Mat) <- c(paste0("X", sprintf("%04d", seq_len(1000))), paste0("P", sprintf("%03d", 1:100)))
  Mat <- data.table(ID = Sigs[, ID], Mat)
  Mat <- merge(Mat, Sigs, by = "ID")
  setkey(Mat, ID)
  Mat <- Mat[Sigs$ID, ]
  fwrite(Mat, paste0("./analysis/21.ABFProcessing/03.MixtureAA/FeatureMatrix/FeatureMatrix_", meta[i, file_name], ".txt"), sep = "\t", quote = F)
})

