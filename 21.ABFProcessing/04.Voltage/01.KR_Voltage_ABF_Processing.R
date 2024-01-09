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


meta <- as.data.table(openxlsx::read.xlsx("./data/Voltage/标准品混合物1.xlsx"))
colnames(meta) <- c("file_name", "amino_acid", "voltage", "concentration", "start_time", "end_time", "baseline_mean")

# Arg(R) 50 mV

File <- file.path("./data/Voltage/21205021.abf")
StartTime <- as.numeric(meta[1, start_time])
EndTime <- as.numeric(meta[1, end_time])

abf <- tryCatch(readABF::readABF(File), error = function(e) NULL)
abf <- as.data.table(as.data.frame(abf))
colnames(abf) <- c("Time", "pA", "mV")[seq_len(ncol(abf))]
abf <- abf[Time > StartTime * 60 & Time < EndTime * 60, ]
abf[, Valid := round(mV) == 50 & pA > 0 & pA < 130]

abf$Sm <- KBins2(sig = abf$pA, pen.value = 1 - 1e-16, minseglen1 = 1, minseglen2 = 1)
saveRDS(abf, "./analysis/21.ABFProcessing/04.Voltage/ABF/ABF_Arg_50mV.Rds")

Sigs <- GetSignal(x = abf, minblockade = 0.1)
Sigs[, DwellTime := EndTime - StartTime]
Sigs <- Sigs[Valid == TRUE]
if(nrow(Sigs) == 0) return(NULL)
Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
Sigs <- data.table(ID = paste(gsub(".abf", "", basename(File)), Sigs[, seq_len(.N)], sep = "_"), Sigs)
fwrite(Sigs, file = paste0("./analysis/21.ABFProcessing/04.Voltage/RawSignal/RawSignal_Arg_50mV.txt"), sep = "\t", row.names = F, quote = F)

Current <- SignalCurrent(Sigs, abf = abf)
saveRDS(Current, file = paste0("./analysis/21.ABFProcessing/04.Voltage/SignalCurrent/SignalCurrent_Arg_50mV.Rds"))



# Arg(R) 75 mV

File <- file.path("./data/Voltage/21205022.abf")
StartTime <- as.numeric(meta[2, start_time])
EndTime <- as.numeric(meta[2, end_time])

abf <- tryCatch(readABF::readABF(File), error = function(e) NULL)
abf <- as.data.table(as.data.frame(abf))
colnames(abf) <- c("Time", "pA", "mV")[seq_len(ncol(abf))]
abf <- abf[Time > StartTime * 60 & Time < EndTime * 60, ]
abf[, Valid := round(mV) == 75 & pA > 0 & pA < 200]

abf$Sm <- KBins2(sig = abf$pA, pen.value = 1 - 1e-16, minseglen1 = 1, minseglen2 = 1)
saveRDS(abf, "./analysis/21.ABFProcessing/04.Voltage/ABF/ABF_Arg_75mV.Rds")

Sigs <- GetSignal(x = abf, minblockade = 0.1)
Sigs[, DwellTime := EndTime - StartTime]
Sigs <- Sigs[Valid == TRUE]
if(nrow(Sigs) == 0) return(NULL)
Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
Sigs <- data.table(ID = paste(gsub(".abf", "", basename(File)), Sigs[, seq_len(.N)], sep = "_"), Sigs)
fwrite(Sigs, file = paste0("./analysis/21.ABFProcessing/04.Voltage/RawSignal/RawSignal_Arg_75mV.txt"), sep = "\t", row.names = F, quote = F)

Current <- SignalCurrent(Sigs, abf = abf)
saveRDS(Current, file = paste0("./analysis/21.ABFProcessing/04.Voltage/SignalCurrent/SignalCurrent_Arg_75mV.Rds"))



# Arg(R) 100 mV

File <- file.path("./data/Voltage/21205022.abf")
StartTime <- as.numeric(meta[3, start_time])
EndTime <- as.numeric(meta[3, end_time])

abf <- tryCatch(readABF::readABF(File), error = function(e) NULL)
abf <- as.data.table(as.data.frame(abf))
colnames(abf) <- c("Time", "pA", "mV")[seq_len(ncol(abf))]
abf <- abf[Time > StartTime * 60 & Time < EndTime * 60, ]
abf[, Valid := round(mV) == 100 & pA > 0 & pA < 280]

abf$Sm <- KBins2(sig = abf$pA, pen.value = 1 - 1e-16, minseglen1 = 1, minseglen2 = 1)
saveRDS(abf, "./analysis/21.ABFProcessing/04.Voltage/ABF/ABF_Arg_100mV.Rds")

Sigs <- GetSignal(x = abf, minblockade = 0.1)
Sigs[, DwellTime := EndTime - StartTime]
Sigs <- Sigs[Valid == TRUE]
if(nrow(Sigs) == 0) return(NULL)
Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
Sigs <- data.table(ID = paste(gsub(".abf", "", basename(File)), Sigs[, seq_len(.N)], sep = "_"), Sigs)
fwrite(Sigs, file = paste0("./analysis/21.ABFProcessing/04.Voltage/RawSignal/RawSignal_Arg_100mV.txt"), sep = "\t", row.names = F, quote = F)

Current <- SignalCurrent(Sigs, abf = abf)
saveRDS(Current, file = paste0("./analysis/21.ABFProcessing/04.Voltage/SignalCurrent/SignalCurrent_Arg_100mV.Rds"))




# Lys(K) 50 mV

File <- file.path("./data/Voltage/21207011.abf")
StartTime <- as.numeric(meta[4, start_time])
EndTime <- as.numeric(meta[4, end_time])

abf <- tryCatch(readABF::readABF(File), error = function(e) NULL)
abf <- as.data.table(as.data.frame(abf))
colnames(abf) <- c("Time", "pA", "mV")[seq_len(ncol(abf))]
abf <- abf[Time > StartTime * 60 & Time < EndTime * 60, ]
abf[, Valid := round(mV) == 50 & pA > 0 & pA < 130]

abf$Sm <- KBins2(sig = abf$pA, pen.value = 1 - 1e-16, minseglen1 = 1, minseglen2 = 1)
saveRDS(abf, "./analysis/21.ABFProcessing/04.Voltage/ABF/ABF_Lys_50mV.Rds")

Sigs <- GetSignal(x = abf, minblockade = 0.1)
Sigs[, DwellTime := EndTime - StartTime]
Sigs <- Sigs[Valid == TRUE]
if(nrow(Sigs) == 0) return(NULL)
Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
Sigs <- data.table(ID = paste(gsub(".abf", "", basename(File)), Sigs[, seq_len(.N)], sep = "_"), Sigs)
fwrite(Sigs, file = paste0("./analysis/21.ABFProcessing/04.Voltage/RawSignal/RawSignal_Lys_50mV.txt"), sep = "\t", row.names = F, quote = F)

Current <- SignalCurrent(Sigs, abf = abf)
saveRDS(Current, file = paste0("./analysis/21.ABFProcessing/04.Voltage/SignalCurrent/SignalCurrent_Lys_50mV.Rds"))



# Lys(K) 75 mV

File <- file.path("./data/Voltage/21207011.abf")
StartTime <- as.numeric(meta[5, start_time])
EndTime <- as.numeric(meta[5, end_time])

abf <- tryCatch(readABF::readABF(File), error = function(e) NULL)
abf <- as.data.table(as.data.frame(abf))
colnames(abf) <- c("Time", "pA", "mV")[seq_len(ncol(abf))]
abf <- abf[Time > StartTime * 60 & Time < EndTime * 60, ]
abf[, Valid := round(mV) == 75 & pA > 0 & pA < 200]

abf$Sm <- KBins2(sig = abf$pA, pen.value = 1 - 1e-16, minseglen1 = 1, minseglen2 = 1)
saveRDS(abf, "./analysis/21.ABFProcessing/04.Voltage/ABF/ABF_Lys_75mV.Rds")

Sigs <- GetSignal(x = abf, minblockade = 0.1)
Sigs[, DwellTime := EndTime - StartTime]
Sigs <- Sigs[Valid == TRUE]
if(nrow(Sigs) == 0) return(NULL)
Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
Sigs <- data.table(ID = paste(gsub(".abf", "", basename(File)), Sigs[, seq_len(.N)], sep = "_"), Sigs)
fwrite(Sigs, file = paste0("./analysis/21.ABFProcessing/04.Voltage/RawSignal/RawSignal_Lys_75mV.txt"), sep = "\t", row.names = F, quote = F)

Current <- SignalCurrent(Sigs, abf = abf)
saveRDS(Current, file = paste0("./analysis/21.ABFProcessing/04.Voltage/SignalCurrent/SignalCurrent_Lys_75mV.Rds"))



# Lys(K) 100 mV

File <- file.path("./data/Voltage/21207011.abf")
StartTime <- as.numeric(meta[6, start_time])
EndTime <- as.numeric(meta[6, end_time])

abf <- tryCatch(readABF::readABF(File), error = function(e) NULL)
abf <- as.data.table(as.data.frame(abf))
colnames(abf) <- c("Time", "pA", "mV")[seq_len(ncol(abf))]
abf <- abf[Time > StartTime * 60 & Time < EndTime * 60, ]
abf[, Valid := round(mV) == 100 & pA > 0 & pA < 280]

abf$Sm <- KBins2(sig = abf$pA, pen.value = 1 - 1e-16, minseglen1 = 1, minseglen2 = 1)
saveRDS(abf, "./analysis/21.ABFProcessing/04.Voltage/ABF/ABF_Lys_100mV.Rds")

Sigs <- GetSignal(x = abf, minblockade = 0.1)
Sigs[, DwellTime := EndTime - StartTime]
Sigs <- Sigs[Valid == TRUE]
if(nrow(Sigs) == 0) return(NULL)
Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
Sigs <- data.table(ID = paste(gsub(".abf", "", basename(File)), Sigs[, seq_len(.N)], sep = "_"), Sigs)
fwrite(Sigs, file = paste0("./analysis/21.ABFProcessing/04.Voltage/RawSignal/RawSignal_Lys_100mV.txt"), sep = "\t", row.names = F, quote = F)

Current <- SignalCurrent(Sigs, abf = abf)
saveRDS(Current, file = paste0("./analysis/21.ABFProcessing/04.Voltage/SignalCurrent/SignalCurrent_Lys_100mV.Rds"))



