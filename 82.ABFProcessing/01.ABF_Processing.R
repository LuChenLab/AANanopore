setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(readABF)
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
  mean(x[x >= min(q) & x <= max(q)])
}

WhiskerRange <- function(x) {
  iqr <- 1.5 * IQR(x, na.rm = TRUE)
  c(min(x[x >= quantile(x, 1/4, na.rm = T) - iqr], na.rm = TRUE), max(x[x <= quantile(x, 3/4, na.rm = T) + iqr], na.rm = TRUE))
}

MainRidge <- function(x, bw = 0.001, bwstep = NULL, n = 1024, peaks = 2, mingap = 0.01, CI = 0.99, plot = T, ...) {
  Peaks <- function(x, y) {
    stopifnot(length(x) == length(y))
    uR <- IRanges(diff(y) >= 0)
    dR <- IRanges(diff(y) <= 0)
    res <- do.call(c, lapply(seq_along(uR), function(i) {
      reduce(c(uR[i], dR[i]))
    }))
    
    Ps <- mapply(seq_along(res), FUN = function(i) {
      x[start(res[i]):end(res[i])][which.max(y[start(res[i]):end(res[i])])]
    })
    Py <- mapply(seq_along(res), FUN = function(i) {
      max(y[start(res[i]):end(res[i])])
    })
    res <- data.table(x = Ps, xmin = x[start(res)], xmax = x[end(res)], y = Py)
    attr(res, "peaks") <- Ps
    return(res)
  }
  den <- density(x, bw = bw, n = n, ...)
  ps <- Peaks(den$x, den$y)
  res <- ps[order(y, decreasing = T)][seq_len(peaks)]
  res <- na.omit(res)
  if(peaks > 1) {
    if(is.null(bwstep)) bwstep <- mingap/10
    while (nrow(res) > 1 & min(res[, diff(sort(x))]) < mingap) {
      bw <- bw + bwstep
      den <- density(x, bw = bw, n = n, ...)
      ps <- Peaks(den$x, den$y)
      res <- ps[order(y, decreasing = T)][seq_len(peaks)]
      res <- na.omit(res)
    }
  }
  if(CI < 1) {
    new_x <- lapply(seq_len(nrow(res)), function(i) {
      resy <- with(den, y[x > res[i, xmin] & x < res[i, xmax]])
      resy <- resy > max(resy) * (1 - CI) + min(resy)
      resx <- range(with(den, x[x > res[i, xmin] & x < res[i, xmax]])[resy])
      resx
    })
    new_x <- do.call(rbind, new_x)
    res$xmin <- new_x[, 1]
    res$xmax <- new_x[, 2]
  }
  res$ratio <- mapply(function(i) mean(x >= res$xmin[i] & x <= res$xmax[i]), seq_len(nrow(res)))
  
  if(plot) {
    plot(den)
    abline(v = res[, x], lty = 2, col = 2)
    abline(v = res[, xmin], lty = 2, col = 3)
    abline(v = res[, xmax], lty = 2, col = 3)
  } else {
    return(res)
  }
}

SignalCurrent <- function(x, abf, cores = 1) {
  target <- abf[inrange(Time, x$StartTime, x$EndTime, incbounds = T)]
  target <- mclapply(seq_len(nrow(x)), function(i) {
    j <- target[between(Time, x[i, StartTime], x[i, EndTime])]
    j[, .(ID = x[i, ID], Current = pA / x[i, BaseMean])]
  }, mc.cores = cores)
  do.call(rbind, target)
}

GetSignal <- function(abf, minblockade = 0.1, cores = 1, minbaseline = 80) {
  abf <- abf[rep(runLength(Rle(abf$Sm)) >= 20, runLength(Rle(abf$Sm)))]
  gr <- IRanges(end = cumsum(runLength(Rle(abf[, Sm]))), width = runLength(Rle(abf[, Sm])))
  
  Sm_N <- abf[Sm > 10, .N, ceiling(Sm)]
  Sm_N <- Sm_N[, .(Sm = ceiling, N, P = N / sum(N) * 100)]
  Sm_N <- Sm_N[P > 1 & Sm > minbaseline]
  Sm_N[, M := Sm * (1 - minblockade)]
  
  S0 <- mclapply(Sm_N[, M], function(x) IRanges(abf[, Sm < x]), mc.cores = cores)
  S0 <- sort(unique(do.call(c, S0)))
  S0 <- S0[start(S0) != 1 & end(S0) != nrow(abf)]
  
  S0 <- IRanges(S0, yL = abf[start(S0) - 1, Sm], yR = abf[end(S0) + 1, Sm], 
                ym = mcmapply(function(x) min(abf[start(S0[x]):end(S0[x]), Sm]), seq_along(S0), mc.cores = cores), 
                yLl = width(gr[subjectHits(findOverlaps(IRanges(start(S0) - 1), gr))]), 
                yRl = width(gr[subjectHits(findOverlaps(IRanges(end(S0) + 1), gr))]), 
                StartTime = abf[start(S0), Time], 
                EndTime = abf[end(S0), Time])
  S0 <- S0[mcols(S0)[, 3] / apply(mcols(S0)[, 1:2], 1, max) < 1 - minblockade] # remove signal blockade too samll
  S0 <- S0[mcols(S0)$yR >= mcols(S0)$yL * (1 - minblockade) & mcols(S0)$yL >= mcols(S0)$yR * (1 - minblockade)] # remove signal with unbalance baseline
  
  Sigs <- mclapply(seq_along(S0), function(j) {
    # print(j)
    x <- c(start(S0[j]) - 1, end(S0[j]) + 1)
    L <- gr[end(gr) == x[1]]
    R <- gr[start(gr) == x[2]]
    DeltaMean <- abs(mean2(abf[start(L):end(L), pA]) - mean2(abf[start(R):end(R), pA]))
    BaselineSD <- tryCatch(sd(c(abf[start(L):end(L), pA], abf[start(R):end(R), pA])), error = function(e) 0)
    if(DeltaMean > BaselineSD) return(NULL) # remove signal with unbalance baseline
    BaseMean <- mean2(c(abf[start(L):end(L), pA], abf[start(R):end(R), pA]))
    if(BaseMean <= minbaseline) return(NULL)
    StartTime <- suppressWarnings(abf[x[1]:x[2]][Sm < BaseMean * (1 - minblockade), min(Time)])
    EndTime <- suppressWarnings(abf[x[1]:x[2]][Sm < BaseMean * (1 - minblockade), max(Time)])
    # StartTime <- with(S0[j], StartTime)
    # EndTime <- with(S0[j], EndTime)
    if(any(is.infinite(c(StartTime, EndTime)))) return(NULL)
    if(abf[between(Time, StartTime, EndTime), any(Sm > BaseMean * (1 - minblockade))]) return(NULL) # remove signal rebound to baseline
    gri <- gr[queryHits(findOverlaps(gr, IRanges(x[1], x[2]), type = "within"))]
    Main_ridge <- tryCatch(abf[between(Time, StartTime, EndTime), MainRidge(x = pA / BaseMean, bw = 0.005, peaks = 1, CI = 0.95, plot = F)], error = function(e) NULL)
    if(is.null(Main_ridge)) return(NULL)
    if(Main_ridge[, x] <= 0) return(NULL)
    StageSD = abf[between(Time, StartTime, EndTime), sd(pA), Sm][, mean(V1)]
    data.table(StartTime = StartTime, EndTime = EndTime, 
               LeftLength = width(L), RightLength = width(R), 
               BaseMean = BaseMean,
               DeltaMean = DeltaMean,
               StageSD = StageSD, 
               CurrentSD = abf[between(Time, StartTime, EndTime), sd(pA)], 
               Segments = max(1, sum(prop.table(width(gri)) > mean(prop.table(width(gri))))), 
               Voltage = abf[between(Time, StartTime, EndTime), mean(mV)], 
               SignalCurrent = BaseMean * Main_ridge[, x], 
               SignalCurrentPercent = Main_ridge$ratio * 100, 
               SignalCurrentWidth = Main_ridge[, xmax - xmin])
  }, mc.cores = cores)
  if(any(mapply(function(x) class(x)[1], Sigs) == "try-error")) return(NULL)
  Sigs <- do.call(rbind, Sigs)
  if(is.null(Sigs)) return(NULL)
  Sigs[, DwellTime := EndTime - StartTime]
  Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
  Sigs <- Sigs[, .SD[which.min(DeltaMean)], .(StartTime, EndTime)]
  return(Sigs)
}

ABF_file <- list.files("./data", "ABF", recursive = T, full.names = T)
ABF_file <- ABF_file[!duplicated(basename(ABF_file))]
names(ABF_file) <- gsub(".ABF", "", basename(ABF_file))

abf_files <- list.files("./data", "abf", recursive = T, full.names = T)
abf_files <- abf_files[!duplicated(basename(abf_files))]
names(abf_files) <- gsub(".abf", "", basename(abf_files))

files <- c(abf_files, ABF_file)
files <- data.table(x = files, y = names(files))
files1 <- do.call(c, lapply(list.files("./data/MetaInfomation", "txt", full.names = T), function(x) fread(x)[[1]]))
files2 <- unlist(lapply(1:8, function(x) as.data.table(openxlsx::read.xlsx("./data/MetaInfomation/整理数据(20230831).xlsx", sheet = x))[[1]]))
files <- files[y %in% c(files1, files2)]
files <- files[!y %in% gsub("FeatureMatrix_", "", gsub(".txt", "", list.files("./analysis/82.ABFProcessing/FeatureMatrix")))]

outdir <- file.path("./analysis/82.ABFProcessing")
files <- files[y %in% gsub(".Rds", "", gsub("ABF_", "", list.files("./analysis/81.ABFProcessing/ABF", ".Rds")))]

StdAA <- fread("./data/MetaInfomation/StandardAA_Meta.txt")
files <- files[y %in% StdAA[[1]]]

lapply(split(seq_len(nrow(files)), cut_number(seq_len(nrow(files)), n = 12))[[12]], function(i) {
  print(paste("V2:", i))
  if(file.exists(file.path(outdir, "FeatureMatrix", paste0("FeatureMatrix_", files[i, y], ".txt")))) return(NULL)
  abf <- readRDS(file.path(gsub("82.ABFProcessing", "81.ABFProcessing", outdir), "ABF", paste0("ABF_", files[i, y], ".Rds")))
  
  # Raw Signal
  if(file.exists(file.path(outdir, "RawSignal", paste0("RawSignal_", files[i, y], ".txt")))) {
    Sigs <- fread(file.path(outdir, "RawSignal", paste0("RawSignal_", files[i, y], ".txt")))
  } else {
    Sigs <- tryCatch(GetSignal(abf = abf, minblockade = 0.1, cores = 10), error = function(e) NULL)
    if(is.null(Sigs)) return(NULL)
    if(nrow(Sigs) == 0) return(NULL)
    Sigs <- data.table(ID = paste(files[i, y], Sigs[, seq_len(.N)], sep = "_"), Sigs)
    Sigs <- Sigs[!is.na(DeltaMean) & !is.na(StageSD) & !is.na(CurrentSD)]
    fwrite(Sigs, file = file.path(outdir, "RawSignal", paste0("RawSignal_", files[i, y], ".txt")), sep = "\t", row.names = F, quote = F)
  }
  
  # Signal Current
  if(file.exists(file.path(outdir, "SignalCurrent", paste0("SignalCurrent_", files[i, y], ".Rds")))) {
    Current <- readRDS(file.path(outdir, "SignalCurrent", paste0("SignalCurrent_", files[i, y], ".Rds")))
  } else {
    Current <- tryCatch(SignalCurrent(Sigs, abf = abf, cores = 10), error = function(e) NULL)
    if(is.null(Current)) return(NULL)
    saveRDS(Current, file = file.path(outdir, "SignalCurrent", paste0("SignalCurrent_", files[i, y], ".Rds")))
  }
  
  # Feature matrix
  Mat <- mclapply(Sigs[, ID], function(j) {
    D <- density(Current[ID == j, Current], from = 0, to = 1, n = 1000, bw = 0.005)$y
    pA <- sort(Current[ID == j, Current])
    pA <- pA[round(seq(1, length(pA), length.out = 100))]
    c(D, pA)
  }, mc.cores = 10)
  Mat <- as.data.table(do.call(rbind, Mat))
  colnames(Mat) <- c(paste0("X", sprintf("%04d", seq_len(1000))), paste0("P", sprintf("%03d", 1:100)))
  Mat <- data.table(ID = Sigs[, ID], Mat)
  Mat <- merge(Mat, Sigs, by = "ID")
  setkey(Mat, ID)
  Mat <- Mat[Sigs$ID, ]
  fwrite(Mat, file = file.path(outdir, "FeatureMatrix", paste0("FeatureMatrix_", files[i, y], ".txt")), sep = "\t", row.names = F, quote = F)
})

