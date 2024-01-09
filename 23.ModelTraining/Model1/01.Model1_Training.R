setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(parallel)
library(Biostrings)

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


files <- list.files("./analysis/22.SignalSelecting/01.StandardAA", full.names = TRUE, recursive = TRUE)
Selected_Sig <- lapply(files, fread)
Selected_Sig <- do.call(rbind, Selected_Sig)
mean(Selected_Sig[, ID] %in% Feature_Mat[, ID])

ggplot(Selected_Sig, aes(x = Blockade, y = DwellTime, label = A, colour = A)) + 
  geom_text() + 
  scale_y_log10()


set.seed(123)
Selected_Sig_Upsample <- Selected_Sig[, .SD[sample(.N, 2000, replace = T), ], A]

ggplot(Selected_Sig_Upsample, aes(x = Blockade, colour = A)) + 
  geom_line(stat = "density")


files <- list.files("./analysis/22.SignalSelecting/02.Background/noise", full.names = TRUE)
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
ggplot(Selected_Sig0_2, aes(x = Blockade, y = DwellTime, label = A, colour = A)) + 
  geom_text() + 
  scale_y_log10()


Selected_Sig_TU <- rbind(Selected_Sig_Upsample, Selected_Sig0_2)
Selected_Sig_TU[, A := plyr::mapvalues(A, Biostrings::AMINO_ACID_CODE, names(Biostrings::AMINO_ACID_CODE))]
Selected_Sig_TU[, A := factor(A, levels = c(AA_STANDARD, "Noise"))]
setkey(Selected_Sig_TU, A)

Train <- setDF(Feature_Mat[Selected_Sig_TU[, ID], 2:1001])
plot(colSums(Train))
Train <- Train[, paste0("X", sprintf("%04d", 401:1000))]
Train$Class <- Selected_Sig_TU[, A]

saveRDS(Train, "./analysis/23.ModelTraining/Model1/Train.Rds")

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

Train <- readRDS("./analysis/23.ModelTraining/Model1/Train.Rds")
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
saveRDS(RFC, file = "./analysis/23.ModelTraining/Model1/RF.Rds")

Train <- subset.data.frame(Train, Class != "C")
Train$Class <- as.character(Train$Class)
levels(Train$Class) <- setdiff(c(Biostrings::AA_STANDARD, "Noise"), "C")

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
saveRDS(RFC, file = "./analysis/23.ModelTraining/Model1/RF_NoCys.Rds")
