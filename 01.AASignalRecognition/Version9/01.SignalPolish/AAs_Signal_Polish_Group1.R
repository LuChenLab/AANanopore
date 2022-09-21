setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(parallel)
library(openxlsx)
library(changepoint)

KBins2 <- function(sig, minseglen1 = 100, minseglen2 = 10, pen.value = 1e-3) {
  ansmean1 <- suppressWarnings(changepoint::cpt.mean(sig, penalty = "MBIC", method = "PELT", minseglen = minseglen1))
  Tab1 <- data.table(P = sig, B1 = rep(seq_len(changepoint::nseg(ansmean1)), changepoint::seg.len(ansmean1)))
  Tab2 <- Tab1[, .(P = median(P)), by = "B1"]
  
  ansmean2 <- suppressWarnings(changepoint::cpt.meanvar(Tab2[, P], penalty = "Asymptotic", pen.value = pen.value, method = "PELT", minseglen = minseglen2))
  Tab2[, B := rep(seq_len(changepoint::nseg(ansmean2)), changepoint::seg.len(ansmean2))]
  res <- merge(Tab1, Tab2[, .(B1, B)], by = "B1")
  res[, B1 := NULL]
  res[, .(P = median(P), N = .N), "B"][, rep(P, N)]
}


meta1 <- data.table(openxlsx::read.xlsx("./data/meta_info_and_base_line_20210525.xlsx"))
meta2 <- data.table(openxlsx::read.xlsx("./data/meta_info_and_base_line_20210525additional.xlsx"))
meta <- rbind(meta1, meta2)
meta[amino_acid == "cys", amino_acid := "Cys"]
meta1 <- meta[!amino_acid %in% c("His", "Lys", "Pro", "Tyr", "Cys", "Gly", "Thr")]
meta1$file_path <- mapply(meta1$file_name, FUN = function(x) list.files("./data", recursive = TRUE, full.names = TRUE, pattern = as.character(x)))

mclapply(seq_len(nrow(meta1)), function(i) {
  if(file.exists(paste0("./analysis/01.AASignalRecognition/Version9/01.SignalPolish/Group1/", meta1[i, file_name], ".Rds"))) return(NULL)
  
  File <- meta1[i, file_path]
  StartTime <- meta1[i, start_time]
  EndTime <- meta1[i, end_time]
  
  abf <- readABF::readABF(File)
  abf <- as.data.table(as.data.frame(abf))
  colnames(abf) <- c("Time", "pA", "mV")
  abf <- abf[Time > StartTime * 60 & Time < EndTime * 60, ]
  
  abf <- abf[round(mV) == 50, ]
  abf <- abf[pA > 0, ]
  abf <- abf[pA < 130, ]
  abf$Sm <- KBins2(sig = abf$pA, pen.value = 1 - 1e-16, minseglen1 = 1, minseglen2 = 1)
  saveRDS(abf, paste0("./analysis/01.AASignalRecognition/Version9/01.SignalPolish/Group1/", meta1[i, file_name], ".Rds"))
}, mc.cores = 2)

