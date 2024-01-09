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

GetSignal <- function(abf, minblockade = 0.1, cores = 1) {
  fun1 <- function(x, y, minblockade = 0.1, cores = cores) {
    y2 <- Rle(y)
    uR <- IRanges(diff(runValue(y2)) >= 0)
    dR <- IRanges(diff(runValue(y2)) <= 0)
    if(start(uR[1]) == 1) uR <- uR[-1]
    if(max(end(dR)) > max(end(uR))) dR <- dR[-c(length(dR))]
    stopifnot(length(uR) == length(dR))
    
    res <- IRanges(end = end(uR), start = start(dR))
    res <- IRanges(start = mcmapply(function(x) sum(runLength(y2)[1:x]) + 1, start(res), mc.cores = cores), 
                   end = mcmapply(function(x) sum(runLength(y2)[1:x]), end(res), mc.cores = cores), 
                   yL = runValue(y2)[start(res)], 
                   yR = runValue(y2)[end(res) + 1],
                   ym = mcmapply(function(x) min(runValue(y2)[start(res[x]):end(res[x])]), seq_along(res), mc.cores = cores), 
                   yLl = runLength(y2)[start(res)], 
                   yRl = runLength(y2)[end(res) + 1])
    res <- res[mcols(res)$ym > 0]
    res <- res[1 - mcols(res)[, 3] / apply(mcols(res)[, 1:2], 1, max) > minblockade] # remove signal blockade too samll
    res <- res[apply(mcols(res)[, 4:5], 1, max) >= 50] # remove signal baseline too short
    ISminblockade <- mcols(res)$yR >= mcols(res)$yL * (1 - minblockade) & mcols(res)$yL >= mcols(res)$yR * (1 - minblockade)
    if(any(!ISminblockade)) {
      res1 <- res[ISminblockade]
      res2 <- res[!ISminblockade]
      mcols(res2)$end2 <- mcmapply(FUN = function(i) {
        ys <- start(IRanges(y2 >= mcols(res2)$yL[i] * (1 - minblockade)))
        min(ys[ys > start(res2[i])]) - 1
      }, seq_along(res2), mc.cores = cores)
      
      mcols(res2)$start2 <- mcmapply(FUN = function(i) {
        ys <- end(IRanges(y2 >= mcols(res2)$yR[i] * (1 - minblockade)))
        max(ys[ys < end(res2[i])]) + 1
      }, seq_along(res2), mc.cores = cores)
      
      res2 <- unique(c(IRanges(start = start(res2), end = ifelse(is.infinite(mcols(res2)$end2), length(x), mcols(res2)$end2)), 
                       IRanges(end = end(res2), start = ifelse(is.infinite(mcols(res2)$start2), 0, mcols(res2)$start2))))
      res2 <- res2[start(res2) != 0 & end(res2) != length(x)]
      res2 <- IRanges(res2, yL = y[start(res2) - 1], yR = y[end(res2) + 1], 
                      ym = mcmapply(function(x) min(y[start(res2[x]):end(res2[x])]), seq_along(res2), mc.cores = cores), 
                      yLl = width(gr[subjectHits(findOverlaps(IRanges(start(res2) - 1), gr))]), 
                      yRl = width(gr[subjectHits(findOverlaps(IRanges(end(res2) + 1), gr))]))
      res2 <- res2[1 - mcols(res2)[, 3] / apply(mcols(res2)[, 1:2], 1, max) > minblockade] # remove signal blockade too samll
      res2 <- res2[mcols(res2)$yR >= mcols(res2)$yL * (1 - minblockade) & mcols(res2)$yL >= mcols(res2)$yR * (1 - minblockade)] # remove signal with unbalance baseline
      res <- sort(c(res1, res2))
    }
    res <- res[apply(mcols(res)[, 4:5], 1, min) >= 50] # remove signal baseline too short
    res <- res[mcols(res)$yR > 50 & mcols(res)$yL > 50] # remove signal baseline too small
    return(res)
  }
  
  gr <- IRanges(end = cumsum(runLength(Rle(abf[, Sm]))), width = runLength(Rle(abf[, Sm])))
  S0 <- abf[, suppressWarnings(fun1(Time, Sm, minblockade = minblockade, cores = cores))]
  
  Sigs <- mclapply(seq_along(S0), function(j) {
    # print(j)
    x <- c(start(S0[j]) - 1, end(S0[j]) + 1)
    L <- gr[end(gr) == x[1]]
    R <- gr[start(gr) == x[2]]
    DeltaMean <- abs(mean2(abf[start(L):end(L), pA]) - mean2(abf[start(R):end(R), pA]))
    BaselineSD <- tryCatch(sd(c(abf[start(L):end(L), pA], abf[start(R):end(R), pA])), error = function(e) 0)
    if(DeltaMean > BaselineSD) return(NULL) # remove signal with unbalance baseline
    BaseMean <- mean2(c(abf[start(L):end(L), pA], abf[start(R):end(R), pA]))
    if(BaseMean <= 0) return(NULL)
    StartTime <- suppressWarnings(abf[x[1]:x[2]][Sm < BaseMean * (1 - minblockade), min(Time)])
    EndTime <- suppressWarnings(abf[x[1]:x[2]][Sm < BaseMean * (1 - minblockade), max(Time)])
    if(any(is.infinite(c(StartTime, EndTime)))) return(NULL)
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
  Sigs[, DwellTime := EndTime - StartTime]
  Sigs[, Blockade := 1 - SignalCurrent / BaseMean]
  return(Sigs)
}

SignalCurrent <- function(x, abf, cores = 1) {
  target <- abf[inrange(Time, x$StartTime, x$EndTime, incbounds = T)]
  target <- mclapply(seq_len(nrow(x)), function(i) {
    j <- target[between(Time, x[i, StartTime], x[i, EndTime])]
    j[, .(ID = x[i, ID], Current = pA / x[i, BaseMean])]
  }, mc.cores = cores)
  do.call(rbind, target)
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
files3 <- rbind(as.data.table(openxlsx::read.xlsx("data/MetaInfomation/neoantigen水解20220314.xlsx", sheet = 1))[, 1:7], 
                as.data.table(openxlsx::read.xlsx("data/ZhangMing_20220406/neoantigen水解20220406.xlsx", sheet = 1))[, 1:7])
files3 <- files3[, unique(file_name)]
files4 <- as.data.table(openxlsx::read.xlsx("data/MetaInfomation/20231020 氨肽酶超滤.xlsx", sheet = 1))[[1]]
files5 <- as.data.table(openxlsx::read.xlsx("data/MetaInfomation/氨基酸浓度梯度实验信息表.xlsx", sheet = 1))[[1]]
files6 <- as.data.table(openxlsx::read.xlsx("data/WangZichun/20231022/20231020起实验数据记录.xlsx", sheet = 1))[[1]]
files7 <- as.data.table(openxlsx::read.xlsx("data/ChenShanchuan/20231022/20231022_Val.xlsx", sheet = 1))[[1]]
files8 <- as.data.table(openxlsx::read.xlsx("data/WangZichun/20231023/20231020起实验数据记录.xlsx", sheet = 1))[[1]]
files9 <- as.data.table(openxlsx::read.xlsx("data/ChenShanchuan/20231024/20231024_data.xlsx", sheet = 1))[[1]]
files10 <- as.data.table(openxlsx::read.xlsx("data/WangZichun/20231025/20231020起实验数据记录.xlsx", sheet = 1))[[1]]
files11 <- as.data.table(openxlsx::read.xlsx("data/20231025实时水解检测/实时水解检测 2.xlsx", sheet = 1))[[1]]
files12 <- as.data.table(openxlsx::read.xlsx("data/ChenShanchuan/20231025/20231025_Arg.xlsx", sheet = 1))[[1]]
files13 <- as.data.table(openxlsx::read.xlsx("data/ChenShanchuan/20231026/20231026_Arg.xlsx", sheet = 1))[[1]]
files14 <- as.data.table(openxlsx::read.xlsx("data/ChenShanchuan/20231026P1/20231026_P1.xlsx", sheet = 1))[[1]]
files15 <- as.data.table(openxlsx::read.xlsx("data/ChenShanchuan/20231026P1/20231026_P1_2.xlsx", sheet = 1))[[1]]
files16 <- as.data.table(openxlsx::read.xlsx("data/WangZichun/20231027/20231020起实验数据记录.xlsx", sheet = 1))[[1]]
files17 <- as.data.table(openxlsx::read.xlsx("data/ChenShanchuan/20231027AcK/20231027_AcK.xlsx", sheet = 1))[[1]]
files18 <- as.data.table(openxlsx::read.xlsx("data/ChenShanchuan/20231027AcK/20231027_AcK2.xlsx", sheet = 1))[[1]]
files19 <- as.data.table(openxlsx::read.xlsx("data/WangZichun/20231028/20231020起实验数据记录.xlsx", sheet = 1))[[1]]
files20 <- as.data.table(openxlsx::read.xlsx("data/WangZichun/20231040CMCC/20231020起实验数据记录.xlsx", sheet = 1))[[1]]
files21 <- as.data.table(openxlsx::read.xlsx("data/ChenShanchuan/20231030/20231030.xlsx", sheet = 1))[[1]]
files22 <- as.data.table(openxlsx::read.xlsx("data/ChenShanchuan/20231031/20231031.xlsx", sheet = 1))[[1]]
files23 <- as.data.table(openxlsx::read.xlsx("data/ChenShanchuan/20231031PTM/20231031ptm.xlsx", sheet = 1))[[1]]
files24 <- as.data.table(openxlsx::read.xlsx("data/ChenShanchuan/20231031PTM2/20231031_PTM2.xlsx", sheet = 1))[[1]]
files25 <- as.data.table(openxlsx::read.xlsx("data/231101实时水解/实时水解检测_实验信息表.xlsx", sheet = 1))[[1]]

files <- files[y %in% c(files1, files2, files3, files4, files5, files6, files7, files8, files9, files10, files11, files12, 
                        files13, files14, files15, files16, files17, files18, files19, files20, files21, files22, files23, files24, files25)]
files <- files[!y %in% gsub("FeatureMatrix_", "", gsub(".txt", "", list.files("./analysis/81.ABFProcessing/FeatureMatrix")))]
files <- files[!y %in% c("21121004", "20231020_0006")]

outdir <- file.path("./analysis/81.ABFProcessing")

lapply(1:nrow(files), function(i) {
  print(i)
  if(file.exists(file.path(outdir, "FeatureMatrix", paste0("FeatureMatrix_", files[i, y], ".txt")))) return(NULL)
  
  if(file.exists(file.path(outdir, "ABF", paste0("ABF_", files[i, y], ".Rds")))) {
    abf <- readRDS(file.path(outdir, "ABF", paste0("ABF_", files[i, y], ".Rds")))
  } else {
    abf <- tryCatch(readABF::readABF(files[i, x]), error = function(e) NULL)
    if(is.null(abf)) return(NULL)
    abf <- as.data.table(as.data.frame(abf))
    colnames(abf) <- c("Time", "pA", "mV")[1:ncol(abf)]
    if(!is.element("mV", colnames(abf))) abf[, mV := NA]
    abf$Sm <- KBins2(sig = abf$pA, pen.value = 1 - 1e-16, minseglen1 = 1, minseglen2 = 1)
    saveRDS(abf, file.path(outdir, "ABF", paste0("ABF_", files[i, y], ".Rds")))
  }
  
  # Raw Signal
  if(file.exists(file.path(outdir, "RawSignal", paste0("RawSignal_", files[i, y], ".txt")))) {
    Sigs <- fread(file.path(outdir, "RawSignal", paste0("RawSignal_", files[i, y], ".txt")))
  } else {
    Sigs <- tryCatch(GetSignal(abf = abf, minblockade = 0.1, cores = 20), error = function(e) NULL)
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
    Current <- tryCatch(SignalCurrent(Sigs, abf = abf, cores = 20), error = function(e) NULL)
    if(is.null(Current)) return(NULL)
    saveRDS(Current, file = file.path(outdir, "SignalCurrent", paste0("SignalCurrent_", files[i, y], ".Rds")))
  }
  
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
  fwrite(Mat, file = file.path(outdir, "FeatureMatrix", paste0("FeatureMatrix_", files[i, y], ".txt")), sep = "\t", row.names = F, quote = F)
})

lapply(split(seq_len(nrow(files)), cut_number(seq_len(nrow(files)), n = 9))[[9]], function(i) {
  print(i)
  if(file.exists(file.path(outdir, "FeatureMatrix", paste0("FeatureMatrix_", files[i, y], ".txt")))) return(NULL)
  
  if(file.exists(file.path(outdir, "ABF", paste0("ABF_", files[i, y], ".Rds")))) {
    abf <- readRDS(file.path(outdir, "ABF", paste0("ABF_", files[i, y], ".Rds")))
  } else {
    abf <- tryCatch(readABF::readABF(files[i, x]), error = function(e) NULL)
    if(is.null(abf)) return(NULL)
    abf <- as.data.table(as.data.frame(abf))
    colnames(abf) <- c("Time", "pA", "mV")[1:ncol(abf)]
    if(!is.element("mV", colnames(abf))) abf[, mV := NA]
    abf$Sm <- KBins2(sig = abf$pA, pen.value = 1 - 1e-16, minseglen1 = 1, minseglen2 = 1)
    saveRDS(abf, file.path(outdir, "ABF", paste0("ABF_", files[i, y], ".Rds")))
  }
  
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