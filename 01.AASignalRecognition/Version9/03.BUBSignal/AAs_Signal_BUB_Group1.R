setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(parallel)
library(openxlsx)
library(S4Vectors)
library(IRanges)

BUB_Sig <- function(Mat, LB = 2000) {
  LRle <- Rle(Mat[, L])
  LRg <- IRanges(start(LRle), end(LRle))
  mcols(LRg)$L <- runValue(LRle)
  
  LRleVa <- paste0(runValue(LRle), collapse = "")
  allcombn <- c("BUB", "BOUB", "BUOB", "BOUOB")
  
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

mean2 <- function(x, q1 = 0.25, q2 = 0.75) {
  q <- quantile(x, c(q1, q2))
  mean(x[x > min(q) & x < max(q)])
}




meta <- do.call(rbind, lapply(list.files("./analysis/01.AASignalRecognition/Version9/02.RangesL0L1/Group1/", ".txt", full.names = TRUE), fread))
meta$minL1 <- apply(meta[, .(MaxL1, MinL1)], 1, min)
meta$maxL1 <- apply(meta[, .(MaxL1, MinL1)], 1, max)
meta[, MinL1 := NULL]
meta[, MaxL1 := NULL]
setnames(meta, "minL1", "MinL1")
setnames(meta, "maxL1", "MaxL1")

mclapply(seq_len(nrow(meta)), function(i) {
  MinL0 <- meta[i, MinL0]
  MaxL0 <- meta[i, MaxL0]
  
  MinL1 <- meta[i, MinL1]
  MaxL1 <- meta[i, MaxL1]
  
  abf <- readRDS(paste0("./analysis/01.AASignalRecognition/Version9/01.SignalPolish/Group1/", meta[i, file_name], ".Rds"))
  
  # BUB
  abf[Sm <= MaxL0 & Sm >= MinL0, L := "B"] # baseline
  abf[Sm < MaxL1, L := "U"] # U shape
  abf[is.na(L), L := "O"] # Other
  
  BUB <- BUB_Sig(Mat = abf)
  BUB <- BUB[mapply(function(x) max(x$Sm), BUB) <= MaxL0]
  
  # Baseline value filtering
  # L0F <- data.table(pA = mapply(function(x) x[L == "B", mean2(pA)], BUB))
  # rg <- quantile(L0F$pA, c(0.05, 0.95))
  # BUB <- BUB[L0F[, pA > min(rg) & pA < max(rg)]]
  
  # Baseline length filtering
  BL_Sum <- mapply(BUB, FUN = function(x) sum(runLength(Rle(x[, L == "B"]))[runValue(Rle(x[, L == "B"])) == TRUE]))
  BL_Min <- mapply(BUB, FUN = function(x) min(runLength(Rle(x[, L == "B"]))[runValue(Rle(x[, L == "B"])) == TRUE]))
  BUB <- BUB[BL_Sum >= 500 & BL_Min >= 50]
  
  # Time filtering
  Wid <- mapply(BUB, FUN = function(x) x[L == "U", diff(range(Time))])
  BUB <- BUB[Wid > 0.00075]
  saveRDS(BUB, paste0("./analysis/01.AASignalRecognition/Version9/03.BUBSignal/Group1/", meta[i, file_name], "_BUB.Rds"))
}, mc.cores = 10)

