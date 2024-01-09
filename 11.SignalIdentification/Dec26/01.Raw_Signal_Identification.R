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
    SignalCurrent <- with(density(abf2[(x[1] + 1):(x[2] - 1), pA], n = 10000), x[which.max(y)])
    SignalCurrentPercent <- mean(abs(abf2[(x[1] + 1):(x[2] - 1), Sm] - SignalCurrent) < 1) * 100
    data.table(StartTime = abf2[x[1], Time], EndTime = abf2[x[2], Time], 
               LeftLength = width(L), RightLength = width(R), 
               BaseMean = mean2(c(abf2[start(L):end(L), pA], abf2[start(R):end(R), pA])),
               DeltaMean = abs(mean2(abf2[start(L):end(L), pA]) - mean2(abf2[start(R):end(R), pA])),
               StageSD = abf2[(x[1] + 1):(x[2] - 1), sd(pA), Sm][, mean(V1)], 
               CurrentSD = abf2[(x[1] + 1):(x[2] - 1), sd(pA)], 
               Segments = sum(prop.table(width(gri)) > 0.1), Valid = all(abf2[start(L):end(R), Valid]), 
               SignalCurrent = SignalCurrent, SignalCurrentPercent = SignalCurrentPercent)
  })
  Sigs <- do.call(rbind, Sigs)
  return(Sigs)
}

SignalCurrent <- function(x, abf, cores = 10) {
  Sigs_Multiple_Current <- mclapply(seq_len(nrow(x)), FUN = function(i) {
    data.table(ID = x[i, ID], Current = abf[Time > x[i, StartTime] & Time < x[i, EndTime], pA] / x[i, BaseMean])
  }, mc.cores = cores)
  do.call(rbind, Sigs_Multiple_Current)
}


meta <- do.call(rbind, lapply(list.files("./analysis/01.AASignalRecognition/Version9/02.RangesL0L1/", ".txt", full.names = TRUE, recursive = T), fread))
meta <- meta[!file_name %in% c("21205012", "21205021", "21201008", "21303003")]
meta <- meta[, .(file_name, amino_acid, concentration, start_time, end_time)]
meta[, file_name := as.character(file_name)]

meta1 <- data.table(openxlsx::read.xlsx("./data/meta_info_and_base_line_20210525.xlsx"))
meta2 <- data.table(openxlsx::read.xlsx("./data/meta_info_and_base_line_20210525additional.xlsx"))
meta <- rbind(meta1, meta2)
meta[amino_acid == "cys", amino_acid := "Cys"]
# meta1 <- meta[!amino_acid %in% c("His", "Lys", "Pro", "Tyr", "Cys", "Gly", "Thr")]
# meta1$file_path <- mapply(meta1$file_name, FUN = function(x) list.files("./data", recursive = TRUE, full.names = TRUE, pattern = as.character(x)))
meta$file_path <- mapply(meta$file_name, FUN = function(x) list.files("./data", recursive = TRUE, full.names = TRUE, pattern = as.character(x)))


mclapply(1:nrow(meta), function(i) {
  print(i)
  File <- meta[i, file_path]
  StartTime <- meta[i, start_time]
  EndTime <- meta[i, end_time]
  
  abf <- readABF::readABF(File)
  abf <- as.data.table(as.data.frame(abf))
  colnames(abf) <- c("Time", "pA", "mV")
  abf <- abf[Time > StartTime * 60 & Time < EndTime * 60, ]
  
  abf[, Valid := round(mV) == 50 & pA > 0 & pA < 130]
  abf$Sm <- KBins2(sig = abf$pA, pen.value = 1 - 1e-16, minseglen1 = 1, minseglen2 = 1)
  Sigs <- GetSignal(x = abf)
  Sigs[, DwellTime := EndTime - StartTime]
  Sigs <- data.table(ID = paste(meta[i, file_name], Sigs[, seq_len(.N)], sep = "_"), Sigs)
  Current <- SignalCurrent(Sigs, abf = abf)
  saveRDS(Current, file = paste0("./analysis/11.SignalIdentification/Dec26/SignalCurrent_", meta[i, file_name], ".Rds"))
  saveRDS(abf, file = paste0("./analysis/11.SignalIdentification/Dec26/ABF_", meta[i, file_name], ".Rds"))
  fwrite(Sigs, file = paste0("./analysis/11.SignalIdentification/Dec26/RawSignal_", meta[i, file_name], ".txt"), sep = "\t", row.names = F, quote = F)
}, mc.cores = 1)


mclapply(1:nrow(meta), function(i) {
  print(i)
  if(file.exists(paste0("./analysis/11.SignalIdentification/Dec26/ABF_", meta[i, file_name], ".Rds")) & file.exists(paste0("./analysis/11.SignalIdentification/Dec26/RawSignal_", meta[i, file_name], ".txt"))) {
    abf <- readRDS(file = paste0("./analysis/11.SignalIdentification/Dec26/ABF_", meta[i, file_name], ".Rds"))
    Sigs <- fread(paste0("./analysis/11.SignalIdentification/Dec26/RawSignal_", meta[i, file_name], ".txt"))
    Sigs[, DwellTime := EndTime - StartTime]
    Sigs <- data.table(ID = paste(meta[i, file_name], Sigs[, seq_len(.N)], sep = "_"), Sigs)
    if(file.exists(paste0("./analysis/11.SignalIdentification/Dec26/SignalCurrent_", meta[i, file_name], ".Rds"))) {
      Current <- readRDS(paste0("./analysis/11.SignalIdentification/Dec26/SignalCurrent_", meta[i, file_name], ".Rds"))
      if(all(Sigs[, ID] %in% Current[, ID])) {
        return(NULL)
      } else {
        Current <- SignalCurrent(Sigs, abf = abf, cores = 1)
        saveRDS(Current, file = paste0("./analysis/11.SignalIdentification/Dec26/SignalCurrent_", meta[i, file_name], ".Rds"))
      }
    } else {
      Current <- SignalCurrent(Sigs, abf = abf, cores = 1)
      saveRDS(Current, file = paste0("./analysis/11.SignalIdentification/Dec26/SignalCurrent_", meta[i, file_name], ".Rds"))
    }
  }
}, mc.cores = 1)



mclapply(1:nrow(meta), function(i) {
  print(i)
  if(file.exists(paste0("./analysis/11.SignalIdentification/Dec26/ABF_", meta[i, file_name], ".Rds"))) {
    abf <- readRDS(file = paste0("./analysis/11.SignalIdentification/Dec26/ABF_", meta[i, file_name], ".Rds"))
    Sigs <- GetSignal(x = abf)
    Sigs[, DwellTime := EndTime - StartTime]
    Sigs <- data.table(ID = paste(meta[i, file_name], Sigs[, seq_len(.N)], sep = "_"), Sigs)
    fwrite(Sigs, file = paste0("./analysis/11.SignalIdentification/Dec26/RawSignal_", meta[i, file_name], "_V2.txt"), sep = "\t", row.names = F, quote = F)
  }
}, mc.cores = 20)



