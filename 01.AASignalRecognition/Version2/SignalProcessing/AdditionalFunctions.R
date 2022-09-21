
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
  Wid <- mapply(BUB, FUN = function(x) x[L == "U", .N]) # 去除时间太长的阻塞以及太短的尖刺
  BUB <- BUB[Wid < quantile(Wid, 0.95) & Wid > quantile(Wid, 0.05) & Wid > 30]
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












