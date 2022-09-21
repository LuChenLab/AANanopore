library(readABF)
library(data.table)
library(S4Vectors)
library(IRanges)
library(changepoint)
library(ggplot2)
library(MASS)
library(strucchange)
library(ggpubr)
library(cowplot)

KBins <- function(sig, minseglen1 = 100, minseglen2 = 10, pen.value = 1e-3) {
  ansmean1 <- suppressWarnings(changepoint::cpt.meanvar(sig, penalty = "None", method = "PELT", minseglen = minseglen1))
  Tab1 <- data.table(P = sig, B1 = rep(seq_len(nseg(ansmean1)), seg.len(ansmean1)))
  Tab2 <- Tab1[, .(P = median(P)), by = "B1"]
  
  ansmean2 <- suppressWarnings(changepoint::cpt.meanvar(Tab2[, P], penalty = "Asymptotic", pen.value = pen.value, method = "PELT", minseglen = minseglen2))
  Tab2[, B := rep(seq_len(nseg(ansmean2)), seg.len(ansmean2))]
  res <- merge(Tab1, Tab2[, .(B1, B)], by = "B1")
  res[, B1 := NULL]
  res[, .(P = median(P), N = .N), "B"][, rep(P, N)]
}

KBins2 <- function(sig, minseglen1 = 100, minseglen2 = 10, pen.value = 1e-3) {
  ansmean1 <- suppressWarnings(changepoint::cpt.mean(sig, penalty = "MBIC", method = "PELT", minseglen = minseglen1))
  Tab1 <- data.table(P = sig, B1 = rep(seq_len(nseg(ansmean1)), seg.len(ansmean1)))
  Tab2 <- Tab1[, .(P = median(P)), by = "B1"]
  
  ansmean2 <- suppressWarnings(changepoint::cpt.meanvar(Tab2[, P], penalty = "Asymptotic", pen.value = pen.value, method = "PELT", minseglen = minseglen2))
  Tab2[, B := rep(seq_len(nseg(ansmean2)), seg.len(ansmean2))]
  res <- merge(Tab1, Tab2[, .(B1, B)], by = "B1")
  res[, B1 := NULL]
  res[, .(P = median(P), N = .N), "B"][, rep(P, N)]
}

KBins3 <- function(sig, minseglen1 = 100, pen.value = 1e-3, penalty = "MBIC") {
  ansmean1 <- suppressWarnings(changepoint::cpt.mean(sig, penalty = penalty, method = "PELT", minseglen = minseglen1, pen.value = pen.value))
  Tab1 <- data.table(P = sig, B1 = rep(seq_len(nseg(ansmean1)), seg.len(ansmean1)))
  Tab2 <- Tab1[, .(Sm = median(P)), by = "B1"]
  res <- merge(Tab1, Tab2, by = "B1")
  res$Sm  
}



L0Coor <- function(den) {
  Bx <- den$x[start(IRanges(diff(den$y) > 0))]
  By <- den$y[start(IRanges(diff(den$y) > 0))]
  
  Px <- den$x[start(IRanges(diff(den$y) < 0))]
  Py <- den$y[start(IRanges(diff(den$y) < 0))]
  
  Py <- Py[Px > 100 & Px < 120]
  Px <- Px[Px > 100 & Px < 120]
  
  Peakx <- which(sort(c(Bx, Px)) == Px[which.max(Py)])
  c(sort(c(Bx, Px))[Peakx - 1], sort(c(Bx, Px))[Peakx + 1])
}

ConfusedL0 <- function(den) {
  Bx <- den$x[start(IRanges(diff(den$y) > 0))]
  By <- den$y[start(IRanges(diff(den$y) > 0))]
  
  Px <- den$x[start(IRanges(diff(den$y) < 0))]
  Py <- den$y[start(IRanges(diff(den$y) < 0))]
  
  Py <- Py[Px > 100 & Px < 120]
  Px <- Px[Px > 100 & Px < 120]
  abs(diff(sort(Py, decreasing = T))[1])/max(Py)
}

BUB_Sig <- function(Mat, LB = 2000) {
  LRle <- Rle(Mat[, L])
  LRg <- IRanges(start(LRle), end(LRle))
  mcols(LRg)$L <- runValue(LRle)
  
  LRleVa <- paste0(runValue(LRle), collapse = "")
  allcombn <- c("BUB", "BOUB", "BUOB")
  
  BUB <- lapply(allcombn, function(x) gregexpr(x, LRleVa)[[1]])
  BUB <- BUB[!mapply(function(x) all(x == -1), BUB)]
  
  BUB <- lapply(BUB, function(x) {
    start <- as.numeric(x)
    steps <- unique(attr(x, "match.length"))
    lapply(start, function(x) {
      rgs <- LRg[x:(x + steps - 1)]
      
      if(width(rgs)[1] > LB) {
        start(rgs)[1] <- end(rgs)[1] - LB
      }
      
      if(width(rgs)[steps] > LB) {
        end(rgs)[steps] <- start(rgs)[steps] + LB
      }
      Mat[start(range(rgs)):end(range(rgs)), ]
    })
  })
  BUB <- do.call(c, BUB)
  
  Wid <- mapply(BUB, FUN = function(x) x[L == "U", diff(range(Time))]) # 去除时间太长的阻塞以及太短的尖刺
  BUB <- BUB[Wid < quantile(Wid, 0.95) & Wid > quantile(Wid, 0.05)]
  BUB
}


L1Coor <- function(den, MinBR = 0.1, MaxBR = 0.3) {
  Bx <- den$x[start(IRanges(diff(den$y) > 0))]
  By <- den$y[start(IRanges(diff(den$y) > 0))]
  
  Px <- den$x[start(IRanges(diff(den$y) < 0))]
  Py <- den$y[start(IRanges(diff(den$y) < 0))]
  
  Py <- Py[Px > MinBR & Px < MaxBR]
  Px <- Px[Px > MinBR & Px < MaxBR]
  
  Peakx <- which(sort(c(Bx, Px)) == Px[which.max(Py)])
  c(sort(c(Bx, Px))[Peakx - 1], sort(c(Bx, Px))[Peakx + 1])
}

library(strucchange)
CISelect <- function(x, adjust = 1, plot = TRUE) {
  densi <- density(x, adjust = adjust)
  denst <- data.table(x = densi$x, y = densi$y)
  fit_bp = breakpoints(y ~ 1, data = denst, breaks = 2)
  cpt <- denst[fit_bp$breakpoints, x]
  pek <- denst[, .SD[which.max(y), x]]
  
  while (all(pek < cpt) | all(pek > cpt)) {
    adjust <- adjust + 0.1
    densi <- density(x, adjust = adjust)
    denst <- data.table(x = densi$x, y = densi$y)
    fit_bp = breakpoints(y ~ 1, data = denst, breaks = 2)
    cpt <- denst[fit_bp$breakpoints, x]
    pek <- denst[, .SD[which.max(y), x]]
  }
  
  if(plot) {
    plot(densi, main = paste0("density(x = x, adjust = ", adjust, ")"))
    abline(v = cpt)
  }
  cpt
}


mean2 <- function(x, q1 = 0.25, q2 = 0.75) {
  q <- quantile(x, c(q1, q2))
  mean(x[x > min(q) & x < max(q)])
}

mad2 <- function(x) {
  median(abs(x - mean(x)))
}

purity <- function(x, w = 2, l = 20) { # here, w is the hight around the main stage, l is the length to include
  sn <- x[L == "U", .N, Sm][N >= l, ]
  mm <- sn[which.max(N), Sm]
  sn[abs(Sm - mm) <= w, Sm := mm]
  max(sn[, sum(N), Sm][, prop.table(V1)])
}

ConfusedL1 <- function(den) {
  Bx <- den$x[start(IRanges(diff(den$y) > 0))]
  By <- den$y[start(IRanges(diff(den$y) > 0))]
  
  Px <- den$x[start(IRanges(diff(den$y) < 0))]
  Py <- den$y[start(IRanges(diff(den$y) < 0))]
  
  # Py <- Py[Px > Min1 & Px < Max1]
  # Px <- Px[Px > Min1 & Px < Max1]
  abs(diff(sort(Py, decreasing = T))[1])/max(Py)
}


plotS <- function(x) {
  ggplot() +
    geom_step(data = x, aes(x = Time, y = pA)) +
    geom_step(data = x, aes(x = Time, y = Sm), color = "red") +
    theme_classic(base_size = 15)
}


plotSig <- function(x) {
  ggplot() + 
    geom_line(data = x, mapping = aes(x = Time, y = pA)) + 
    geom_line(data = x[L == "U"], mapping = aes(x = Time, y = pA), colour = "red", size = 1.1) + 
    theme_minimal(base_size = 15)
}

plotSig2 <- function(x) {
  ggplot() + 
    geom_line(data = x, mapping = aes(x = Time, y = pA)) + 
    # geom_line(data = x[L == "U"], mapping = aes(x = Time, y = pA), colour = "red", size = 1.1) + 
    geom_line(data = x, mapping = aes(x = Time, y = Sm), colour = "red", size = 1.1) + 
    theme_minimal(base_size = 15)
}


mean2 <- function(x, q1 = 0.25, q2 = 0.75, step = 0.05) {
  q <- quantile(x, c(q1, q2))
  while (sum(x > min(q) & x < max(q)) < 1) {
    q1 <- q1 - step
    q2 <- q2 + step
    q <- quantile(x, c(q1, q2))
  }
  mean(x[x > min(q) & x < max(q)])
}

mean3 <- function(x, HT = 10) {
  if(length(x) >= 5*HT) {
    mean(x[(HT+1):(length(x) - HT)])
  } else {
    HT <- floor(length(x)/5)
    mean(x[(HT+1):(length(x) - HT)])
  }
}

BinSignal <- function(x) {
  Sigs <- x[L == "U", unique(Sm)]
  do.call(rbind, lapply(Sigs, function(s) {
    DwellTime <- x[L == "U" & Sm == s, diff(range(Time)) * 1000]
    StartTime <- x[L == "U" & Sm == s, min(Time)]
    BaseLineMean <- x[L == "B", mean2(pA)]
    SignalMean <- x[L == "U" & Sm == s, mean3(pA)]
    SignalSD <- x[L == "U" & Sm == s, sd(pA)]
    Blockade <- 1 - SignalMean/BaseLineMean
    data.table(StartTime = StartTime, DwellTime, BaseLineMean, SignalMean, SignalSD, Blockade)
  }))
}

L0Coor <- function(den, d = NULL) {
  Bx <- den$x[start(IRanges(diff(den$y) > 0))]
  By <- den$y[start(IRanges(diff(den$y) > 0))]
  
  Px <- den$x[start(IRanges(diff(den$y) < 0))]
  Py <- den$y[start(IRanges(diff(den$y) < 0))]
  
  Py <- Py[Px > 100 & Px < 120]
  Px <- Px[Px > 100 & Px < 120]
  
  if(is.null(d)) {
    Peakx <- which(sort(c(Bx, Px)) == Px[which.max(Py)])
    c(sort(c(Bx, Px))[Peakx - 1], sort(c(Bx, Px))[Peakx + 1])
  } else {
    Peakx <- which(sort(c(Bx, Px)) == Px[which.min(abs(Px - d))])
    c(sort(c(Bx, Px))[Peakx - 1], sort(c(Bx, Px))[Peakx + 1])
  }
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

DensityPeaksSignal <- function(Sig, adjust = 0.5, plot = TRUE) {
  BaseLineMean <- Sig[L == "B", mean2(pA)]
  L1 <- BaseLineMean * c(1 - c(0.11, 0.26))
  
  den <- density(x = Sig[pA < max(L1), pA], adjust = adjust)
  Pks <- tryCatch(PeaksOfDensity(den), error = function(e) NULL)
  while (is.null(Pks)) {
    adjust <- adjust + 0.1
    den <- density(x = Sig[pA < max(L1), pA], adjust = adjust)
    Pks <- tryCatch(PeaksOfDensity(den), error = function(e) NA)
  }
  
  if(!is.null(L1)) {
    Pks <- Pks[Px > min(L1) & Px < max(L1)]
  }
  
  while (nrow(Pks) > 1 & Pks[, max(Py)] / Pks[, sort(Py, decreasing = TRUE)][2] < 10) {
    adjust <- adjust + 0.1
    den <- density(x = Sig[pA < max(L1), pA], adjust = adjust)
    Pks <- PeaksOfDensity(den)
    
    if(!is.null(L1)) {
      Pks <- Pks[Px > min(L1) & Px < max(L1)]
    }
  }
  
  Pks <- Pks[which.max(Py)]
  
  StartTime <- Sig[L == "U", min(Time)]
  TimeRatio <- Sig[L == "U", mean(Sm > Pks[, B1x] & Sm < Pks[, B2x])]
  DwellTime <- Sig[L == "U", diff(range(Time)) * 1000] * TimeRatio
  AllTime <- Sig[L == "U", diff(range(Time)) * 1000]
  
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
    data.table(StartTime, AllTime, TimeRatio, DwellTime, BaseLineMean, SignalMean, SignalSD, Blockade)
  }
}
