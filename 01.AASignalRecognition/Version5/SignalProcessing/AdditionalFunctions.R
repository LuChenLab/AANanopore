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

DensityCrossArea <- function(x1, x2, adjust = 1, step = 0.1, plot = TRUE, col.x = 1, col.y = 2, lty.x = 1, lty.y = 2) {
  if(min(x1) > max(x2) | max(x1) < min(x2)) {
    # No any overlaps
    d1 <- density(x1, adjust = adjust)
    d2 <- density(x2, adjust = adjust)
    
    areas <- 0
    
    if(plot) {
      xl <- range(c(d1$x, d2$x))
      yl <- c(0, max(c(d1$y, d2$y)))
      plot(d1, xlim = xl, ylim = yl, col = col.x, lty = lty.x,
           main = paste0("N = ", length(x), ", ", length(y), "; adjust = ", adjust),
           xlab = NA)
      lines(d2, xlim = xl, ylim = yl, col = col.y, lty = lty.y)
    }
  }
  
  if(all(x1 >= min(x2) & x1 <= max(x2)) | all(x2 >= min(x1) & x2 <= max(x1))) {
    # Someone within another
    
    if(all(x1 >= min(x2) & x1 <= max(x2))) {
      # x1 within x2
      
      Fn1 <- ecdf(x1) # a *function*
      Fn2 <- ecdf(x2) # a *function*
      
      # Crossover point1
      ad <- adjust
      d1 <- density(x1, adjust = ad)
      d2 <- density(x2, adjust = ad)
      
      fun1 <- function(x) density(x1, from = x, to = x, n = 1, adjust = ad)$y
      fun2 <- function(x) density(x2, from = x, to = x, n = 1, adjust = ad)$y
      solution1 <- tryCatch(uniroot(function(x) fun1(x) - fun2(x),
                                    lower = min(c(x1, x2)),
                                    upper = d1$x[which.max(d1$y)],
                                    tol = 1e-100,
                                    maxiter = 10000)$root, error = function(e) NA)
      
      while (is.na(solution1)) {
        ad <- ad + step
        d1 <- density(x1, adjust = ad)
        d2 <- density(x2, adjust = ad)
        
        fun1 <- function(x) density(x1, from = x, to = x, n = 1, adjust = ad)$y
        fun2 <- function(x) density(x2, from = x, to = x, n = 1, adjust = ad)$y
        solution1 <- tryCatch(uniroot(function(x) fun1(x) - fun2(x),
                                      lower = min(c(x1, x2)),
                                      upper = d1$x[which.max(d1$y)],
                                      tol = 1e-100,
                                      maxiter = 10000)$root, error = function(e) NA)
      }
      
      # Crossover point2
      solution2 <- tryCatch(uniroot(function(x) fun1(x) - fun2(x),
                                    lower = d1$x[which.max(d1$y)],
                                    upper = max(c(x1, x2)),
                                    tol = 1e-100,
                                    maxiter = 10000)$root, error = function(e) NA)
      while (is.na(solution2)) {
        ad <- ad + step
        d1 <- density(x1, adjust = ad)
        d2 <- density(x2, adjust = ad)
        
        fun1 <- function(x) density(x1, from = x, to = x, n = 1, adjust = ad)$y
        fun2 <- function(x) density(x2, from = x, to = x, n = 1, adjust = ad)$y
        solution2 <- tryCatch(uniroot(function(x) fun1(x) - fun2(x),
                                      lower = d1$x[which.max(d1$y)],
                                      upper = max(c(x1, x2)),
                                      tol = 1e-100,
                                      maxiter = 10000)$root, error = function(e) NA)
      }
      
      areas <- Fn1(solution1) +  (1 - Fn1(solution2))
      
      if(plot) {
        xl <- range(c(d1$x, d2$x))
        yl <- c(0, max(c(d1$y, d2$y)))
        
        plot(d1, xlim = xl, ylim = yl, col = col.x, lty = lty.x,
             main = paste0("N = ", length(x), ", ", length(y), "; adjust = ", ad),
             xlab = NA)
        lines(d2, xlim = xl, ylim = yl, col = col.y, lty = lty.y)
        
        polygon(c(solution2, d1$x[d1$x >= solution2]), c(solution2, d1$y[d1$x >= solution2]),
                col = rgb(1, 0, 0, alpha = 0.5),
                border = "red")
        
        polygon(c(d1$x[d1$x <= solution1], solution1), c(d1$y[d1$x <= solution1], solution1),
                col = rgb(1, 0, 0, alpha = 0.5),
                border = "red")
      }
      
    } else {
      # x2 within x1
      Fn1 <- ecdf(x1) # a *function*
      Fn2 <- ecdf(x2) # a *function*
      
      # Crossover point1
      ad <- adjust
      d1 <- density(x1, adjust = ad)
      d2 <- density(x2, adjust = ad)
      
      fun1 <- function(x) density(x1, from = x, to = x, n = 1, adjust = ad)$y
      fun2 <- function(x) density(x2, from = x, to = x, n = 1, adjust = ad)$y
      solution1 <- tryCatch(uniroot(function(x) fun1(x) - fun2(x),
                                    lower = min(c(x1, x2)),
                                    upper = d2$x[which.max(d2$y)],
                                    tol = 1e-100,
                                    maxiter = 10000)$root, error = function(e) NA)
      
      while (is.na(solution1)) {
        ad <- ad + step
        d1 <- density(x1, adjust = ad)
        d2 <- density(x2, adjust = ad)
        
        fun1 <- function(x) density(x1, from = x, to = x, n = 1, adjust = ad)$y
        fun2 <- function(x) density(x2, from = x, to = x, n = 1, adjust = ad)$y
        solution1 <- tryCatch(uniroot(function(x) fun1(x) - fun2(x),
                                      lower = min(c(x1, x2)),
                                      upper = d2$x[which.max(d2$y)],
                                      tol = 1e-100,
                                      maxiter = 10000)$root, error = function(e) NA)
      }
      
      # Crossover point2
      solution2 <- tryCatch(uniroot(function(x) fun1(x) - fun2(x),
                                    lower = d2$x[which.max(d2$y)],
                                    upper = max(c(x1, x2)),
                                    tol = 1e-100,
                                    maxiter = 10000)$root, error = function(e) NA)
      while (is.na(solution2)) {
        ad <- ad + step
        d1 <- density(x1, adjust = ad)
        d2 <- density(x2, adjust = ad)
        
        fun1 <- function(x) density(x1, from = x, to = x, n = 1, adjust = ad)$y
        fun2 <- function(x) density(x2, from = x, to = x, n = 1, adjust = ad)$y
        solution2 <- tryCatch(uniroot(function(x) fun1(x) - fun2(x),
                                      lower = d2$x[which.max(d2$y)],
                                      upper = max(c(x1, x2)),
                                      tol = 1e-100,
                                      maxiter = 10000)$root, error = function(e) NA)
      }
      
      areas <- Fn1(solution2) - Fn1(solution1)
      
      if(plot) {
        xl <- range(c(d1$x, d2$x))
        yl <- c(0, max(c(d1$y, d2$y)))
        
        plot(d1, xlim = xl, ylim = yl, col = col.x, lty = lty.x,
             main = paste0("N = ", length(x), ", ", length(y), "; adjust = ", ad),
             xlab = NA)
        lines(d2, xlim = xl, ylim = yl, col = col.y, lty = lty.y)
        
        polygon(c(solution1, d1$x[d1$x >= solution1 & d1$x <= solution2], solution2), c(0, d1$y[d1$x >= solution1 & d1$x <= solution2], 0),
                col = rgb(1, 0, 0, alpha = 0.5),
                border = "red")
      }
    }
  }
  
  if((min(x1) > min(x2) & max(x1) > max(x2)) | (min(x2) > min(x1) & max(x2) > max(x1))) {
    # Overlap
    
    Fn1 <- ecdf(x1) # a *function*
    Fn2 <- ecdf(x2) # a *function*
    
    ad <- adjust
    d1 <- density(x1, adjust = ad)
    d2 <- density(x2, adjust = ad)
    fun1 <- function(x) density(x1, from = x, to = x, n = 1, adjust = ad)$y
    fun2 <- function(x) density(x2, from = x, to = x, n = 1, adjust = ad)$y
    solution <- tryCatch(uniroot(function(x) fun1(x) - fun2(x), interval = range(c(x1, x2)), tol = 1e-100, maxiter = 10000)$root, error = function(e) NA)
    
    while (is.na(solution)) {
      ad <- ad + step
      d1 <- density(x1, adjust = ad)
      d2 <- density(x2, adjust = ad)
      fun1 <- function(x) density(x1, from = x, to = x, n = 1, adjust = ad)$y
      fun2 <- function(x) density(x2, from = x, to = x, n = 1, adjust = ad)$y
      solution <- tryCatch(uniroot(function(x) fun1(x) - fun2(x), interval = range(c(x1, x2)), tol = 1e-100, maxiter = 10000)$root, error = function(e) NA)
    }
    areas <- min(Fn1(solution), 1 - Fn1(solution))
    
    if(plot) {
      xl <- range(c(d1$x, d2$x))
      yl <- c(0, max(c(d1$y, d2$y)))
      
      plot(d1, xlim = xl, ylim = yl, col = col.x, lty = lty.x,
           main = paste0("N = ", length(x), ", ", length(y), "; adjust = ", ad),
           xlab = NA)
      lines(d2, xlim = xl, ylim = yl, col = col.y, lty = lty.y)
      
      if(Fn1(solution) < 1 - Fn1(solution)) {
        polygon(c(d1$x[d1$x <= solution], solution), c(d1$y[d1$x <= solution], solution),
                col = rgb(1, 0, 0, alpha = 0.5),
                border = "red")
      } else {
        polygon(c(solution, d1$x[d1$x >= solution]), c(solution, d1$y[d1$x >= solution]),
                col = rgb(1, 0, 0, alpha = 0.5),
                border = "red")
      }
    }
  }
  return(areas)
}


PeakCoor <- function(x, adjust = 1) {
  denres <- density(x, adjust = adjust)
  data.table(x = denres$x[which.max(denres$y)], y = max(denres$y))
}
