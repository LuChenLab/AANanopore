setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(ggplot2)
library(parallel)
library(caret)
library(openxlsx)

# Polypeptide1

Signal1 <- as.data.table(read.xlsx("./analysis/02.PolypeptideSequencing/20211025/Version2/Polypeptide1_APRLRFYSL.xlsx"))
Signal1$ID <- paste0(Signal1[, Sample], "_", sprintf("%04d", do.call(c, lapply(Signal1[, .N, Sample][, N], function(x) seq_len(x)))))

Models <- list.files("./analysis/03.MachineLearning/Version7/02.Modeling/01.RF/01.Models", "Rds", full.names = T)

Mat <- readRDS("./analysis/03.MachineLearning/Version7/03.PredictPolypeptide/01.PolypeptideSignalTransformation/Polypeptide1_APRLRFYSL.signal.Rds")
Polypeptide1 <- mclapply(1:6, function(i) {
  Fit <- readRDS(Models[i])
  if(nrow(subset.data.frame(Mat, AllTime > i & SignalSD < 4)) == 0) return(NULL)
  ROC_Test <- predict(Fit, subset.data.frame(Mat, AllTime > i & SignalSD < 4), type = "prob")
  ppm <- melt(as.data.table(ROC_Test, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
  ppm[, Model := paste0("RF", i)]
  ppm <- merge(ppm, as.data.table(Mat, keep.rownames = "ID")[, .(ID, AllTime, DwellTime, SignalSD, Blockade)], by = "ID")
  return(ppm)
}, mc.cores = 6)
Polypeptide1 <- do.call(rbind, Polypeptide1)
Polypeptide1[, pred := as.character(pred)]
Polypeptide1[, AA := plyr::mapvalues(pred, AMINO_ACID_CODE, names(AMINO_ACID_CODE))]

Polypeptide1 <- merge(Polypeptide1, Signal1[, .(ID, TimeRatio)], by = "ID")

ggplot(Polypeptide1[Model == "RF3" & Prob > 0.7 & TimeRatio > 0.8 & DwellTime > 0.75 & SignalSD < 2], aes(x = Blockade)) + 
  geom_histogram(binwidth = 0.002) + 
  theme_bw(base_size = 16)

Polypeptide1[Model == "RF3" & Prob > 0.7 & TimeRatio > 0.8 & DwellTime > 0.75 & SignalSD < 2, .N, AA]
Polypeptide1[Model == "RF6" & Prob > 0.7, .N, AA]

ggplot(Polypeptide1[Prob > 0.8], aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_text(aes(label = AA)) +
  facet_wrap( ~ Model) + 
  scale_color_manual(values = c(RColorBrewer::brewer.pal(n = 12, "Paired"), RColorBrewer::brewer.pal(n = 8, "Dark2")))

Polypeptide1[Model == "RF6" & Prob > 0.5, .N, AA][order(N, decreasing = T)]




# Polypeptide2

Signal2 <- as.data.table(read.xlsx("./analysis/02.PolypeptideSequencing/20211025/Version2/Polypeptide2_RPVKVYPNGAEDESAEAFPLEF.xlsx"))
Signal2$ID <- paste0(Signal2[, Sample], "_", sprintf("%04d", do.call(c, lapply(Signal2[, .N, Sample][, N], function(x) seq_len(x)))))

Mat <- readRDS("./analysis/03.MachineLearning/Version7/03.PredictPolypeptide/01.PolypeptideSignalTransformation/Polypeptide2_RPVKVYPNGAEDESAEAFPLEF.signal.Rds")
Polypeptide2 <- mclapply(1:6, function(i) {
  Fit <- readRDS(Models[i])
  if(nrow(subset.data.frame(Mat, AllTime > i & SignalSD < 4)) == 0) return(NULL)
  ROC_Test <- predict(Fit, subset.data.frame(Mat, AllTime > i & SignalSD < 4), type = "prob")
  ppm <- melt(as.data.table(ROC_Test, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
  ppm[, Model := paste0("RF", i)]
  ppm <- merge(ppm, as.data.table(Mat, keep.rownames = "ID")[, .(ID, AllTime, DwellTime, SignalSD, Blockade)], by = "ID")
  return(ppm)
}, mc.cores = 6)
Polypeptide2 <- do.call(rbind, Polypeptide2)
Polypeptide2[, pred := as.character(pred)]
Polypeptide2[, AA := plyr::mapvalues(pred, AMINO_ACID_CODE, names(AMINO_ACID_CODE))]

Polypeptide2 <- merge(Polypeptide2, Signal2[, .(ID, TimeRatio)], by = "ID")
Polypeptide2[Model == "RF2" & Prob > 0.75, .N, AA][order(N, decreasing = T)]

ggplot(Polypeptide2[Model == "RF3" & Prob > 0.7 & TimeRatio > 0.8 & DwellTime > 0.75 & SignalSD < 2], aes(x = Blockade)) + 
  geom_histogram(binwidth = 0.002) + 
  theme_bw(base_size = 16)

ggplot(Polypeptide2[Prob > 0.8], aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_text(aes(label = AA)) +
  facet_wrap( ~ Model) + 
  scale_color_manual(values = c(RColorBrewer::brewer.pal(n = 12, "Paired"), RColorBrewer::brewer.pal(n = 8, "Dark2")))

Polypeptide2[Model == "RF3" & Prob > 0.5, .N, AA][order(N, decreasing = T)]








# Polypeptide3

Signal3 <- as.data.table(read.xlsx("./analysis/02.PolypeptideSequencing/20211025/Version2/Polypeptide3_DRVYIHPFHL.xlsx"))
Signal3$ID <- paste0(Signal3[, Sample], "_", sprintf("%04d", do.call(c, lapply(Signal3[, .N, Sample][, N], function(x) seq_len(x)))))

Mat <- readRDS("./analysis/03.MachineLearning/Version7/03.PredictPolypeptide/01.PolypeptideSignalTransformation/Polypeptide3_DRVYIHPFHL.signal.Rds")
Polypeptide3 <- mclapply(1:6, function(i) {
  Fit <- readRDS(Models[i])
  if(nrow(subset.data.frame(Mat, AllTime > i & SignalSD < 4)) == 0) return(NULL)
  ROC_Test <- predict(Fit, subset.data.frame(Mat, AllTime > i & SignalSD < 4), type = "prob")
  ppm <- melt(as.data.table(ROC_Test, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
  ppm[, Model := paste0("RF", i)]
  ppm <- merge(ppm, as.data.table(Mat, keep.rownames = "ID")[, .(ID, AllTime, DwellTime, SignalSD, Blockade)], by = "ID")
  return(ppm)
}, mc.cores = 6)
Polypeptide3 <- do.call(rbind, Polypeptide3)
Polypeptide3[, pred := as.character(pred)]
Polypeptide3[, AA := plyr::mapvalues(pred, AMINO_ACID_CODE, names(AMINO_ACID_CODE))]

Polypeptide3 <- merge(Polypeptide3, Signal3[, .(ID, TimeRatio)], by = "ID")
Polypeptide3[Model == "RF3" & Prob > 0.75, .N, AA][order(N, decreasing = T)]


ggplot(Polypeptide3[Prob > 0.8], aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_text(aes(label = AA)) +
  facet_wrap( ~ Model) + 
  scale_color_manual(values = c(RColorBrewer::brewer.pal(n = 12, "Paired"), RColorBrewer::brewer.pal(n = 8, "Dark2")))

Polypeptide3[Model == "RF3" & Prob > 0.5, .N, AA][order(N, decreasing = T)]


