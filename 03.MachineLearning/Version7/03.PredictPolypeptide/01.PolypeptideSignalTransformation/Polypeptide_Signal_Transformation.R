# Polypeptide
setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(openxlsx)

mean2 <- function(x, q1 = 0.25, q2 = 0.75) {
  q <- quantile(x, c(q1, q2))
  mean(x[x >= min(q) & x <= max(q)])
}

# Polypeptide1 APRLRFYSL

Signal <- as.data.table(read.xlsx("./analysis/02.PolypeptideSequencing/20211025/Version2/Polypeptide1_APRLRFYSL.xlsx"))
Signal$ID <- paste0(Signal[, Sample], "_", sprintf("%04d", do.call(c, lapply(Signal[, .N, Sample][, N], function(x) seq_len(x)))))

BUBs <- paste0("./analysis/02.PolypeptideSequencing/20211025/Version2/BUB", Signal[, unique(Sample)], ".Rds")
BUBs <- lapply(BUBs, readRDS)
BUBs <- do.call(c, BUBs)
names(BUBs) <- Signal$ID

BUBs <- mclapply(BUBs, function(x) {
  basemean <- x[L == "B", mean2(pA)]
  x[, pA := pA / basemean]
  x[, Sm := NULL]
  x[L == "U"]
}, mc.cores = 10)

Signal$AllTime <- mapply(BUBs, FUN = function(x) x[, diff(range(Time))] * 1000)

BinExp <- mclapply(BUBs, function(x) round(density(x[, pA], from = 0, to = 1, n = 200, adjust = 0.5)$y, 3), mc.cores = 10)
BinExp <- do.call(rbind, BinExp)
colnames(BinExp) <- paste0("X", sprintf("%03d", seq_len(ncol(BinExp))))
Mat <- merge(as.data.table(BinExp, keep.rownames = "ID"), Signal[, .(ID, AllTime, DwellTime, SignalSD, Blockade)], by = "ID")
Mat <- data.frame(Mat[, -1], row.names = Mat[[1]])

saveRDS(Mat, "./analysis/03.MachineLearning/Version7/03.PredictPolypeptide/01.PolypeptideSignalTransformation/Polypeptide1_APRLRFYSL.signal.Rds")




# Polypeptide2 RPVKVYPNGAEDESAEAFPLEF

Signal <- as.data.table(read.xlsx("./analysis/02.PolypeptideSequencing/20211025/Version2/Polypeptide2_RPVKVYPNGAEDESAEAFPLEF.xlsx"))
Signal$ID <- paste0(Signal[, Sample], "_", sprintf("%04d", do.call(c, lapply(Signal[, .N, Sample][, N], function(x) seq_len(x)))))

BUBs <- paste0("./analysis/02.PolypeptideSequencing/20211025/Version2/BUB", Signal[, unique(Sample)], ".Rds")
BUBs <- lapply(BUBs, readRDS)
BUBs <- do.call(c, BUBs)
names(BUBs) <- Signal$ID

BUBs <- mclapply(BUBs, function(x) {
  basemean <- x[L == "B", mean2(pA)]
  x[, pA := pA / basemean]
  x[, Sm := NULL]
  x[L == "U"]
}, mc.cores = 10)

Signal$AllTime <- mapply(BUBs, FUN = function(x) x[, diff(range(Time))] * 1000)

BinExp <- mclapply(BUBs, function(x) round(density(x[, pA], from = 0, to = 1, n = 200, adjust = 0.5)$y, 3), mc.cores = 10)
BinExp <- do.call(rbind, BinExp)
colnames(BinExp) <- paste0("X", sprintf("%03d", seq_len(ncol(BinExp))))
Mat <- merge(as.data.table(BinExp, keep.rownames = "ID"), Signal[, .(ID, AllTime, DwellTime, SignalSD, Blockade)], by = "ID")
Mat <- data.frame(Mat[, -1], row.names = Mat[[1]])

saveRDS(Mat, "./analysis/03.MachineLearning/Version7/03.PredictPolypeptide/01.PolypeptideSignalTransformation/Polypeptide2_RPVKVYPNGAEDESAEAFPLEF.signal.Rds")



# Polypeptide3 DRVYIHPFHL

Signal <- as.data.table(read.xlsx("./analysis/02.PolypeptideSequencing/20211025/Version2/Polypeptide3_DRVYIHPFHL.xlsx"))
Signal$ID <- paste0(Signal[, Sample], "_", sprintf("%04d", do.call(c, lapply(Signal[, .N, Sample][, N], function(x) seq_len(x)))))

BUBs <- paste0("./analysis/02.PolypeptideSequencing/20211025/Version2/BUB", Signal[, unique(Sample)], ".Rds")
BUBs <- lapply(BUBs, readRDS)
BUBs <- do.call(c, BUBs)
names(BUBs) <- Signal$ID

BUBs <- mclapply(BUBs, function(x) {
  basemean <- x[L == "B", mean2(pA)]
  x[, pA := pA / basemean]
  x[, Sm := NULL]
  x[L == "U"]
}, mc.cores = 10)

Signal$AllTime <- mapply(BUBs, FUN = function(x) x[, diff(range(Time))] * 1000)

BinExp <- mclapply(BUBs, function(x) round(density(x[, pA], from = 0, to = 1, n = 200, adjust = 0.5)$y, 3), mc.cores = 10)
BinExp <- do.call(rbind, BinExp)
colnames(BinExp) <- paste0("X", sprintf("%03d", seq_len(ncol(BinExp))))
Mat <- merge(as.data.table(BinExp, keep.rownames = "ID"), Signal[, .(ID, AllTime, DwellTime, SignalSD, Blockade)], by = "ID")
Mat <- data.frame(Mat[, -1], row.names = Mat[[1]])

saveRDS(Mat, "./analysis/03.MachineLearning/Version7/03.PredictPolypeptide/01.PolypeptideSignalTransformation/Polypeptide3_DRVYIHPFHL.signal.Rds")



