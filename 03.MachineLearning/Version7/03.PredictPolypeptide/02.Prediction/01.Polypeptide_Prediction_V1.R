setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(ggplot2)
library(parallel)
library(caret)
library(openxlsx)

Models <- list.files("./analysis/03.MachineLearning/Version7/02.Modeling/01.RF/01.Models", "Rds", full.names = T)
RF1 <- readRDS(Models[1])
RF2 <- readRDS(Models[2])
RF3 <- readRDS(Models[3])
RF4 <- readRDS(Models[4])
RF5 <- readRDS(Models[5])
RF6 <- readRDS(Models[6])

# Polypeptide1

Signal1 <- as.data.table(read.xlsx("./analysis/02.PolypeptideSequencing/20211025/Version2/Polypeptide1_APRLRFYSL.xlsx"))
Signal1$ID <- paste0(Signal1[, Sample], "_", sprintf("%04d", do.call(c, lapply(Signal1[, .N, Sample][, N], function(x) seq_len(x)))))
Signal1 <- Signal1[TimeRatio > 0.8 & DwellTime > 0.75 & SignalSD < 2]

Mat <- readRDS("./analysis/03.MachineLearning/Version7/03.PredictPolypeptide/01.PolypeptideSignalTransformation/Polypeptide1_APRLRFYSL.signal.Rds")
Mat <- Mat[row.names(Mat) %in% Signal1$ID, ]

Polypeptide1 <- mclapply(1:6, function(i) {
  Fit <- readRDS(Models[i])
  ROC_Test <- predict(Fit, Mat, type = "prob")
  ppm <- melt(as.data.table(ROC_Test, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
  ppm[, Model := paste0("RF", i)]
  ppm <- merge(ppm, as.data.table(Mat, keep.rownames = "ID")[, .(ID, AllTime, DwellTime, SignalSD, Blockade)], by = "ID")
  return(ppm)
}, mc.cores = 6)

Polypeptide1 <- do.call(rbind, Polypeptide1)
Polypeptide1[, pred := as.character(pred)]
Polypeptide1[, AA := plyr::mapvalues(pred, AMINO_ACID_CODE, names(AMINO_ACID_CODE))]

Polypeptide1[Model == "RF6" & Prob > 0.7, .N, AA][order(N, decreasing = T)]

ggplot(Polypeptide1, aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_point() + 
  facet_wrap( ~ Model) + 
  scale_color_manual(values = c(RColorBrewer::brewer.pal(n = 12, "Paired"), RColorBrewer::brewer.pal(n = 8, "Dark2")))


ggplot(Polypeptide1, aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_text(aes(label = AA)) +
  facet_wrap( ~ Model) + 
  scale_color_manual(values = c(RColorBrewer::brewer.pal(n = 12, "Paired"), RColorBrewer::brewer.pal(n = 8, "Dark2")))





# Polypeptide2

Signal2 <- as.data.table(read.xlsx("./analysis/02.PolypeptideSequencing/20211025/Version2/Polypeptide2_RPVKVYPNGAEDESAEAFPLEF.xlsx"))
Signal2$ID <- paste0(Signal2[, Sample], "_", sprintf("%04d", do.call(c, lapply(Signal2[, .N, Sample][, N], function(x) seq_len(x)))))
Signal2 <- Signal2[TimeRatio > 0.8 & SignalSD < 2 & DwellTime > 0.75]

Mat <- readRDS("./analysis/03.MachineLearning/Version7/03.PredictPolypeptide/01.PolypeptideSignalTransformation/Polypeptide2_RPVKVYPNGAEDESAEAFPLEF.signal.Rds")
Mat <- Mat[row.names(Mat) %in% Signal2$ID, ]

Polypeptide2 <- mclapply(1:6, function(i) {
  Fit <- readRDS(Models[i])
  ROC_Test <- predict(Fit, Mat, type = "prob")
  ppm <- melt(as.data.table(ROC_Test, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
  ppm[, Model := paste0("RF", i)]
  ppm <- merge(ppm, as.data.table(Mat, keep.rownames = "ID")[, .(ID, AllTime, DwellTime, SignalSD, Blockade)], by = "ID")
  return(ppm)
}, mc.cores = 6)

Polypeptide2 <- do.call(rbind, Polypeptide2)
Polypeptide2[, pred := as.character(pred)]
Polypeptide2[, AA := plyr::mapvalues(pred, AMINO_ACID_CODE, names(AMINO_ACID_CODE))]

Polypeptide2[Model == "RF2" & Prob > 0.7, .N, AA][order(N, decreasing = T)]

ggplot(Polypeptide2, aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_point() + 
  facet_wrap( ~ Model) + 
  scale_color_manual(values = c(RColorBrewer::brewer.pal(n = 12, "Paired"), RColorBrewer::brewer.pal(n = 8, "Dark2")))






# Polypeptide3

Signal3 <- as.data.table(read.xlsx("./analysis/02.PolypeptideSequencing/20211025/Version2/Polypeptide3_DRVYIHPFHL.xlsx"))
Signal3$ID <- paste0(Signal3[, Sample], "_", sprintf("%04d", do.call(c, lapply(Signal3[, .N, Sample][, N], function(x) seq_len(x)))))
Signal3 <- Signal3[TimeRatio > 0.8 & SignalSD < 2 & DwellTime > 0.75]

Mat <- readRDS("./analysis/03.MachineLearning/Version7/03.PredictPolypeptide/01.PolypeptideSignalTransformation/Polypeptide3_DRVYIHPFHL.signal.Rds")
Mat <- Mat[row.names(Mat) %in% Signal3$ID, ]

Polypeptide3 <- mclapply(1:6, function(i) {
  Fit <- readRDS(Models[i])
  ROC_Test <- predict(Fit, Mat, type = "prob")
  ppm <- melt(as.data.table(ROC_Test, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
  ppm[, Model := paste0("RF", i)]
  ppm <- merge(ppm, as.data.table(Mat, keep.rownames = "ID")[, .(ID, AllTime, DwellTime, SignalSD, Blockade)], by = "ID")
  return(ppm)
}, mc.cores = 6)

Polypeptide3 <- do.call(rbind, Polypeptide3)
Polypeptide3[, pred := as.character(pred)]
Polypeptide3[, AA := plyr::mapvalues(pred, AMINO_ACID_CODE, names(AMINO_ACID_CODE))]

Polypeptide3[Model == "RF2" & Prob > 0.7, .N, AA][order(N, decreasing = T)]

ggplot(Polypeptide3, aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_point() + 
  facet_wrap( ~ Model) + 
  scale_color_manual(values = c(RColorBrewer::brewer.pal(n = 12, "Paired"), RColorBrewer::brewer.pal(n = 8, "Dark2"))) + 
  scale_y_sqrt()

