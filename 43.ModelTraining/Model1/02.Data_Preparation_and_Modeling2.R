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



RFC <- readRDS(file = "./analysis/43.ModelTraining/Model1/RF.Rds")

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
cfM <- cfM/colSums(cfM)
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
  geom_text(aes(label = AA, colour = AA)) + 
  scale_y_log10()



Sigs999 <- rbind(Sigs[Group == "GSA" & D > .9], 
                 Sigs[(Group == "TNRK" & D > .95 & AA != "K") | (Group == "TNRK" & D > .95 & AA == "K" & Blockade > 0.17) ], 
                 Sigs[Group == "QVML" & D > .95], 
                 Sigs[Group == "IYDPFW" & D > .95], 
                 Sigs[Group == "EH" & D > .9])

ggplot(Sigs999, aes(x = Blockade, y = DwellTime)) + 
  geom_text(aes(label = AA, colour = AA)) + 
  scale_y_log10()


merge(Sigs[Pred != "Noise" & Prob > .9 & Pred == AA, .N, AA], Sigs999[, .N, AA], by = "AA")

Sigs[AA == "Y", .N]
Sigs[Pred == "Y" & Pred == AA, ]

ggplot(Sigs[Pred %in% c("I", "Y", "D") & Pred == AA, ], aes(x = Blockade, y = DwellTime)) + 
  geom_text(aes(label = AA, colour = AA)) + 
  scale_y_log10()

ggplot(Sigs[Pred %in% c("I", "Y", "D") & Pred == AA & Prob > .4, ], aes(x = Blockade, colour = AA)) + 
  geom_density()


ggplot(Sigs[Pred %in% c("N", "R", "K") & Pred == AA & Prob > .8, ], aes(x = Blockade, y = DwellTime)) + 
  geom_text(aes(label = AA, colour = AA)) + 
  scale_y_log10()

ggplot(Sigs[Pred %in% c("N", "R", "K") & Pred == AA & Prob > .8, ], aes(x = Blockade, colour = AA)) + 
  geom_density()




Train2 <- unique(rbind(Sigs[Pred != "Noise" & AA == Pred & Prob > .9, .(ID, AA)], Sigs999[, .(ID, AA)]))
saveRDS(Train2, "./analysis/43.ModelTraining/Model1/Train_Data_ID2.Rds")

set.seed(123)
Selected_Sig_Upsample <- Train2[, .SD[sample(.N, 2000, replace = T), ], AA]


files <- list.files("./analysis/22.SignalSelecting/02.Background/noise", full.names = TRUE)
Selected_Sig0 <- lapply(files, fread)
Selected_Sig0 <- do.call(rbind, Selected_Sig0)
Selected_Sig0[, AA := "Noise"]
Selected_Sig0 <- Selected_Sig0[ID %in% Feature_Mat0$ID]




Selected_Sig_TU <- rbind(Selected_Sig_Upsample, Selected_Sig0[, .(ID, AA)], fill = TRUE)
Selected_Sig_TU[, AA := factor(AA, levels = c(setdiff(AA_STANDARD, "C"), "Noise"))]
setkey(Selected_Sig_TU, AA)

Train <- setDF(Feature_Mat[Selected_Sig_TU[, ID], ])
Train <- Train[, c(paste0("X", sprintf("%04d", 601:1000)), "StageSD", "CurrentSD", "Segments", "SignalCurrentPercent", "SignalCurrentWidth", "Blockade", "DwellTime")]
Train$Class <- Selected_Sig_TU[, AA]

saveRDS(Train, "./analysis/43.ModelTraining/Model1/Train2.Rds")

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

Train <- readRDS("./analysis/43.ModelTraining/Model1/Train2.Rds")
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
saveRDS(RFC, file = "./analysis/43.ModelTraining/Model1/RF2.Rds")






