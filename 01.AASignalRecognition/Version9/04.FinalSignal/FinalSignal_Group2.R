setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(parallel)
library(openxlsx)
library(S4Vectors)
library(IRanges)

mean2 <- function(x, q1 = 0.25, q2 = 0.75) {
  q <- quantile(x, c(q1, q2))
  mean(x[x >= min(q) & x <= max(q)])
}

PeaksOfDensity <- function(den) {
  Bx <- den$x[start(IRanges(diff(den$y) > 0))]
  By <- den$y[start(IRanges(diff(den$y) > 0))]
  
  Px <- den$x[start(IRanges(diff(den$y) < 0))]
  Py <- den$y[start(IRanges(diff(den$y) < 0))]
  
  if(min(Px) == min(den$x)) {
    Px <- Px[-1]
    Py <- Py[-1]
  }
  
  if(max(Px) > max(Bx)) {
    Bx <- c(Bx, tail(den$x, 1))
    By <- c(By, tail(den$y, 1))
  }
  
  if(min(Px) < min(Bx)) {
    Bx <- c(head(den$x, 1), Bx)
    By <- c(head(den$y, 1), By)
  }
  
  B1x <- mapply(Px, FUN = function(i) max(Bx[Bx < i]))
  B1y <- By[mapply(B1x, FUN = function(i) which(Bx == i))]
  
  B2x <- mapply(Px, FUN = function(i) min(Bx[Bx > i]))
  B2y <- By[mapply(B2x, FUN = function(i) which(Bx == i))]
  
  Pks <- data.table(Px = Px, Py = Py, B1x = B1x, B1y = B1y, B2x = B2x, B2y = B2y)
}

DensityPeaksSignal <- function(Sig, MinL1, MaxL1, L2min, L2max, adjust = 0.5, plot = TRUE) {
  BaseLineMean <- Sig[L == "B", mean2(pA)]
  L1 <- c(MinL1, MaxL1)
  L2 <- BaseLineMean * (1 - c(L2max, L2min))
  
  Fn1 <- ecdf(Sig[L == "U", pA])
  AreaRatio_L1 <- 1 - Fn1(MinL1)/Fn1(MaxL1)
  AreaRatio_L2 <- ifelse(Fn1(max(L2)) == 0, 0, 1 - Fn1(min(L2))/Fn1(max(L2)))
  
  den <- density(x = Sig[L == "U", pA], adjust = adjust)
  Pks <- tryCatch(PeaksOfDensity(den), error = function(e) NULL)
  while (is.null(Pks)) {
    adjust <- adjust + 0.1
    den <- density(x = Sig[L == "U", pA], adjust = adjust)
    Pks <- tryCatch(PeaksOfDensity(den), error = function(e) NA)
  }
  
  if(!is.null(L1)) {
    Pks <- Pks[Px > min(L1) & Px < max(L1)]
  }
  
  while (nrow(Pks) > 1 & Pks[, max(Py)] / Pks[, sort(Py, decreasing = TRUE)][2] < 10) {
    adjust <- adjust + 0.1
    den <- density(x = Sig[L == "U", pA], adjust = adjust)
    Pks <- PeaksOfDensity(den)
    
    if(!is.null(L1)) {
      Pks <- Pks[Px > min(L1) & Px < max(L1)]
    }
  }
  
  Pks <- Pks[which.max(Py)]
  
  StartTime <- Sig[L == "U", min(Time)]
  TimeRatio <- Sig[L == "U", mean(Sm > Pks[, B1x] & Sm < Pks[, B2x])]
  DwellTime <- Sig[L == "U", diff(range(Time)) * 1000] * TimeRatio
  
  if(Sig[, sum(L == "U" & pA > Pks[, B1x] & pA < Pks[, B2x])] >= 20) {
    SignalMean <- Sig[L == "U" & pA > Pks[, B1x] & pA < Pks[, B2x], mean2(pA)]
    SignalSD <- Sig[L == "U" & pA > Pks[, B1x] & pA < Pks[, B2x], sd(pA)]
    
    if(Sig[, sum(L == "U" & pA > SignalMean)] < 20) {
      SignalMean <- NA
      SignalSD <- NA
    }
  } else {
    SignalMean <- NA
    SignalSD <- NA
  }
  
  Blockade <- 1 - SignalMean/BaseLineMean
  Outer <- (SignalMean - Sig[, min(Sm)])/(BaseLineMean - SignalMean)
  
  if(plot) {
    par(mfrow = c(2, 2))
    Sig[, plot(Time, pA, type = "s")]
    Sig[, lines(Time, Sm, type = "s", col = 2, lwd = 1.5)]
    abline(h = L1, col = 3)
    
    plot(den, main = paste0("Density, adjust = ", adjust))
    Pks[, points(x = Px, Py, col = 2, pch = 16)]
    Pks[, points(x = B1x, B1y, col = 3, pch = 16)]
    Pks[, points(x = B2x, B2y, col = 3, pch = 16)]
    abline(v = L1, col = 3)
    
    Sig[, plot(Time, pA, type = "s")]
    Sig[, lines(Time, Sm, type = "s", col = 2, lwd = 1.5)]
    abline(h = Pks[, B1x, B2x], col = 3)
    
    Sig[, plot(Time, pA, type = "s", main = paste0("Blockade = ", round(Blockade, 4), ", DwellTime = ", round(DwellTime, 3)))]
    Sig[, lines(Time, Sm, type = "s", col = 2, lwd = 1.5)]
    abline(h = SignalMean, col = 2)
  } else {
    data.table(StartTime, TimeRatio, AreaRatio_L1 = AreaRatio_L1, AreaRatio_L2 = AreaRatio_L2, DwellTime, BaseLineMean, SignalMean, SignalSD, Blockade, Outer = Outer)
  }
}

DensityPeaksSignal <- function(Sig, MinL1, MaxL1, L2min, L2max, adjust = 0.5, plot = TRUE) {
  BaseLineMean <- Sig[L == "B", mean2(pA)]
  L1 <- c(MinL1, MaxL1)
  L2 <- BaseLineMean * (1 - c(L2max, L2min))
  
  Fn1 <- ecdf(Sig[L == "U", Sm])
  AreaRatio_L1 <- Fn1(MaxL1) - Fn1(MinL1)
  AreaRatio_L2 <- ifelse(Fn1(max(L2)) == 0, 0, Fn1(max(L2)) - Fn1(min(L2)))
  Outer <- Fn1(min(L2))
  
  StartTime <- Sig[L == "U", min(Time)]
  DwellTime <- Sig[L == "U", diff(range(Time)) * 1000] * AreaRatio_L1
  
  SignalMean <- Sig[L == "U" & Sm > min(L1) & Sm < max(L1), mean2(pA)]
  SignalSD <- Sig[L == "U" & Sm > min(L1) & Sm < max(L1), sd(pA)]
  
  Blockade <- 1 - SignalMean/BaseLineMean
  data.table(StartTime, AreaRatio_L1, AreaRatio_L2, Outer, DwellTime, BaseLineMean, SignalMean, SignalSD, Blockade)
}


meta <- do.call(rbind, lapply(list.files("./analysis/01.AASignalRecognition/Version9/02.RangesL0L1/Group2/", ".txt", full.names = TRUE), fread))
meta$minL1 <- apply(meta[, .(MaxL1, MinL1)], 1, min)
meta$maxL1 <- apply(meta[, .(MaxL1, MinL1)], 1, max)
meta[, MinL1 := NULL]
meta[, MaxL1 := NULL]
setnames(meta, "minL1", "MinL1")
setnames(meta, "maxL1", "MaxL1")

aa_L2 <- as.data.table(read.xlsx("./data/aa_L2.xlsx"))
meta <- merge(meta, aa_L2, by.x = "amino_acid", by.y = "name")

mclapply(seq_len(nrow(meta)), function(i) {
  print(i)
  BUB <- readRDS(paste0("./analysis/01.AASignalRecognition/Version9/03.BUBSignal/Group2/", meta[i, file_name], "_BUB.Rds"))
  
  Sigs <- lapply(BUB, FUN = function(x) DensityPeaksSignal(Sig = x, MinL1 = meta[i, MinL1], MaxL1 = meta[i, MaxL1], L2min = meta[i, L2min], L2max = meta[i, L2max]))
  Sigs <- do.call(rbind, Sigs)
  fwrite(Sigs, paste0("./analysis/01.AASignalRecognition/Version9/04.FinalSignal/Group2/", meta[i, file_name], "_Sigs.txt"), sep = "\t", row.names = F, quote = F)
}, mc.cores = 1)

