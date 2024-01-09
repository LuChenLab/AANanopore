setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(parallel)
library(Biostrings)

files2 <- list.files("./analysis/13.MachineLearning/01.DataPreparation/FeatureMatrix/Density200Points/BackgroundSignal", pattern = ".txt", full.names = TRUE)
Feature_Mat0 <- lapply(files2, fread)
Feature_Mat0 <- do.call(rbind, Feature_Mat0)
Feature_Mat0 <- na.omit(Feature_Mat0)

files <- list.files("./analysis/13.MachineLearning/01.DataPreparation/FeatureMatrix/Density200Points/StandardAA", pattern = ".txt", full.names = TRUE)
Feature_Mat <- lapply(files, fread)
Feature_Mat <- do.call(rbind, Feature_Mat)

Feature_Mat <- rbind(Feature_Mat, Feature_Mat0)
setkey(Feature_Mat, ID)


files <- list.files("./analysis/11.SignalIdentification/Jan07/AASignal", pattern = "Selected_V2", full.names = TRUE)
Selected_Sig <- lapply(files, fread)
Selected_Sig <- do.call(rbind, Selected_Sig)
mean(Selected_Sig[, ID] %in% Feature_Mat[, ID])
Selected_Sig <- Selected_Sig[A != "Cys"]

files <- list.files("./analysis/12.SignalFiltering/AASignal/Shiny3/Cys", "RawSignal_", full.names = T)
CysSignals <- do.call(rbind, lapply(files, fread))
Cys_ID <- fread("./analysis/12.SignalFiltering/AASignal/Shiny3/Cys/Cys_ID.txt")
CysSignals <- CysSignals[ID %in% Cys_ID[N > 1, ID]]

Selected_Sig <- rbind(Selected_Sig, CysSignals)
Selected_Sig[, .N, A]

ggplot(Selected_Sig, aes(x = Blockade, y = DwellTime, label = A, colour = A)) + 
  geom_text() + 
  scale_y_log10()

set.seed(123)
Selected_Sig_Upsample <- Selected_Sig[, .SD[sample(.N, 2000, replace = T), ], A]

ggplot(Selected_Sig_Upsample, aes(x = Blockade, y = DwellTime, label = A, colour = A)) + 
  geom_text() + 
  scale_y_log10()


files <- list.files("./analysis/11.SignalIdentification/Jan07/BackgroundSignal", pattern = "Selected", full.names = TRUE)
Selected_Sig0 <- lapply(files, fread)
Selected_Sig0 <- do.call(rbind, Selected_Sig0)
Selected_Sig0[, A := "Noise"]
Selected_Sig0 <- Selected_Sig0[ID %in% Feature_Mat0$ID]


ggplot(Selected_Sig0, aes(x = Blockade, y = DwellTime, label = A, colour = A)) + 
  geom_text() + 
  scale_y_log10()

Selected_Sig0[, Blockade2 := round(Blockade * 100)]
set.seed(123)
Selected_Sig0_2 <- rbind(Selected_Sig0[, .SD[sample(.N, 100, replace = T), ], Blockade2], 
                         Selected_Sig0[Blockade < 0.3, .SD[sample(.N, 200, replace = T), ], Blockade2])
Selected_Sig0_2[, Blockade2 := NULL]
ggplot(Selected_Sig0_2, aes(x = Blockade)) + 
  geom_histogram()


Selected_Sig_TU <- rbind(Selected_Sig_Upsample, Selected_Sig0_2)
Selected_Sig_TU[, A := plyr::mapvalues(A, Biostrings::AMINO_ACID_CODE, names(Biostrings::AMINO_ACID_CODE))]
Selected_Sig_TU[, A := factor(A, levels = c(AA_STANDARD, "Noise"))]
setkey(Selected_Sig_TU, A)

Train <- setDF(Feature_Mat[Selected_Sig_TU[, ID], 1:300])
plot(colSums(Train))
Train <- Train[, paste0("X", sprintf("%03d", 81:200))]
Train$Class <- Selected_Sig_TU[, A]

saveRDS(Train, "./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model4/Train.Rds")

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
cl <- makePSOCKcluster(10)
registerDoParallel(cl)

Train <- readRDS("./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model4/Train.Rds")
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
saveRDS(RFC, file = "./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model4/RF_Fit.Rds")

