setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(ggplot2)
library(parallel)
library(caret)

get_density <- function(x, y, n = 1000, ...) {
  dens <- MASS::kde2d(x, y, n = n, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii] / max(dens$z[ii]))
}

AAGroup <- list(GSA = c('G', 'S', 'A'), 
                TNRK = c("T", "N", "R", "K"), 
                QVML = c("Q", "V", "M", "L"), 
                IYDPFW = c("I", "Y", "D", "P", "F", "W"), 
                EH = c("E", "H"))
AAGroup <- data.table(Group = rep(names(AAGroup), mapply(length, AAGroup)), AA = unlist(AAGroup))
AAGroup[, amino_acid := plyr::mapvalues(AA, names(AMINO_ACID_CODE), AMINO_ACID_CODE)]

Sigs <- lapply(AAGroup[, amino_acid], function(aat) {
  sigs <- do.call(rbind, lapply(list.files(paste0("./analysis/42.SignalSelecting/01.StandardAA/", aat), full.names = T), fread))
  sigs$File <- stringr::str_remove_all(sigs$ID, "_([[:digit:]]+)$")
  sigs$D <- sigs[, get_density(x = Blockade, y = log10(DwellTime))]
  sigs
})
Sigs <- do.call(rbind, Sigs)
Sigs <- merge(Sigs, AAGroup, by.x = "A", by.y = "amino_acid")
Sigs[, AA := plyr::mapvalues(A, AMINO_ACID_CODE, names(AMINO_ACID_CODE))]
Sigs[, AA := factor(AA, levels = c('G', 'S', 'A', 'T', 'N', 'R', 'K', 'V', 'Q', 'L', 'M', 'I', 'Y', 'D', 'P', 'F', 'W', 'E', 'H', 'Noise'))]



RFC <- readRDS(file = "./analysis/43.ModelTraining/Model1/RF3.Rds")

files <- list.files("./analysis/21.ABFProcessing/02.Background/FeatureMatrix", pattern = "", full.names = TRUE, recursive = TRUE)
Feature_Mat0 <- lapply(files, fread)
Feature_Mat0 <- do.call(rbind, Feature_Mat0)
Feature_Mat0 <- na.omit(Feature_Mat0)

files <- list.files("./analysis/21.ABFProcessing/01.StandardAA/FeatureMatrix", pattern = "", full.names = TRUE, recursive = TRUE)
Feature_Mat <- lapply(files, fread)
Feature_Mat <- do.call(rbind, Feature_Mat)

setkey(Feature_Mat0, ID)
Feature_Mat <- rbind(Feature_Mat, Feature_Mat0)
setkey(Feature_Mat, ID)

Sigs$Pred <- predict(RFC, Feature_Mat[Sigs$ID, ])
Sigs$Prob <- apply(predict(RFC, Feature_Mat[Sigs$ID, ], type = "prob"), 1, max)
Sigs[, Pred := factor(Pred, levels = c('G', 'S', 'A', 'T', 'N', 'R', 'K', 'V', 'Q', 'L', 'M', 'I', 'Y', 'D', 'P', 'F', 'W', 'E', 'H', 'Noise'))]

cfM <- Sigs[Prob > .9, confusionMatrix(AA, Pred)]
cfM <- as.matrix(cfM)[1:19, 1:19]
cfM <- t(t(cfM)/colSums(cfM))
pheatmap::pheatmap(cfM, cluster_rows = F, cluster_cols = F)



ggplot(Sigs, aes(x = AA, y = Prob)) + 
  geom_violin() + 
  geom_boxplot(width = .1)

ggplot(Sigs, aes(x = Pred, y = Prob)) + 
  geom_violin()

Sigs[, .N, Pred]
Sigs[Pred != "Noise", .N, Pred]
Sigs[Pred != "Noise" & Prob > .9, .N, Pred]
Sigs[Pred != "Noise" & Prob > .9, mean(AA == Pred)]



ggplot(Sigs[Pred != "Noise" & Prob > .8 & Pred == AA, ], aes(x = Blockade, y = DwellTime)) + 
  geom_text(aes(label = AA, colour = AA), size = 2) + 
  scale_y_log10()

ggplot(Sigs[Pred != "Noise" & Pred == AA, ], aes(x = Blockade, y = DwellTime)) + 
  geom_text(aes(label = AA, colour = AA), size = 2) + 
  scale_y_log10()


ggplot(Sigs[D > .8], aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_point() + 
  scale_y_log10()

ggplot(Sigs[Pred != "Noise" & Pred == AA, ], aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_point(alpha = .4) + 
  scale_y_log10() + 
  scale_colour_manual(values = AA_Cols)
#

Train4 <- Sigs[Pred != "Noise", .(ID, AA)]

saveRDS(Train4, "./analysis/43.ModelTraining/Model1/Train_Data_ID4.Rds")

set.seed(123)
Selected_Sig_Upsample <- rbind(Train4[AA %in% Train4[, .N, AA][N < 1000, AA], .SD[sample(.N, 1200, replace = T), ], AA], 
                               Train4[AA %in% Train4[, .N, AA][N >= 1000, AA], ])
setkey(Selected_Sig_Upsample, AA)
Selected_Sig_Upsample[, AA := factor(AA, levels = c('G', 'S', 'A', 'T', 'N', 'R', 'K', 'V', 'Q', 'L', 'M', 'I', 'Y', 'D', 'P', 'F', 'W', 'E', 'H', 'Noise'))]


files <- list.files("./analysis/22.SignalSelecting/02.Background/noise", full.names = TRUE)
Selected_Sig0 <- lapply(files, fread)
Selected_Sig0 <- do.call(rbind, Selected_Sig0)
Selected_Sig0[, AA := "Noise"]
Selected_Sig0 <- Selected_Sig0[ID %in% Feature_Mat0$ID]




Selected_Sig_TU <- rbind(Selected_Sig_Upsample, Selected_Sig0[, .(ID, AA)], fill = TRUE)
Selected_Sig_TU[, AA := factor(AA, levels = c('G', 'S', 'A', 'T', 'N', 'R', 'K', 'V', 'Q', 'L', 'M', 'I', 'Y', 'D', 'P', 'F', 'W', 'E', 'H', 'Noise'))]
setkey(Selected_Sig_TU, AA)




Train <- setDF(Feature_Mat[Selected_Sig_TU[, ID], ])
Train <- Train[, c(paste0("X", sprintf("%04d", 601:1000)), "StageSD", "CurrentSD", "Segments", "SignalCurrentPercent", "SignalCurrentWidth", "Blockade", "DwellTime")]
Train$Class <- Selected_Sig_TU[, AA]

saveRDS(Train, "./analysis/43.ModelTraining/Model1/Train4.Rds")

setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(ggplot2)
library(parallel)
library(caret)

fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10)

library(doParallel)
cl <- makePSOCKcluster(20)
registerDoParallel(cl)

Train <- readRDS("./analysis/43.ModelTraining/Model1/Train4.Rds")
set.seed(825)
RFC <- train(Class ~ ., data = Train, 
             # preProc = c("center", "scale", "YeoJohnson", "nzv"), 
             method = "rf", 
             trControl = fitControl,
             verbose = FALSE,
             ## to evaluate:
             tuneGrid = expand.grid(mtry = 30),
             # tuneLength = 50,
             metric = "Accuracy", 
             allowParallel = TRUE)
saveRDS(RFC, file = "./analysis/43.ModelTraining/Model1/RF4.Rds")

