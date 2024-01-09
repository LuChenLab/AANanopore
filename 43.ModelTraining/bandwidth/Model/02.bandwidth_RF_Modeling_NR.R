setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(ggplot2)
library(parallel)
library(caret)
load("./analysis/43.ModelTraining/bandwidth/TrainTestSet/TrainTestSet_NR.RData")

fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10)

library(doParallel)
cl <- makePSOCKcluster(20)
registerDoParallel(cl)

lapply(c(1e-05, 1:9/10000, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1), function(d) {
  print(d)
  FeaMat <- do.call(rbind, lapply(list.files(paste0("./analysis/43.ModelTraining/bandwidth/FeatureMatrix/bandwidth", d*1000), full.names = T), fread))
  setkey(FeaMat, ID)
  Train <- setDF(FeaMat[Train_Set[, ID]])
  Train <- Train[, c(paste0("X", sprintf("%04d", 601:1000)), "StageSD", "CurrentSD", "Segments", "SignalCurrentPercent", "SignalCurrentWidth", "Blockade", "DwellTime")]
  Train$Class <- factor(Train_Set$AA, levels = c('N', 'R'))
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
  saveRDS(RFC, file = paste0("./analysis/43.ModelTraining/bandwidth/Model/NR_RF", d*1000, ".Rds"))
})


Res <- lapply(c(1e-05, 1:9/10000, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1), function(d) {
  print(d)
  FeaMat <- do.call(rbind, lapply(list.files(paste0("./analysis/43.ModelTraining/bandwidth/FeatureMatrix/bandwidth", d*1000), full.names = T), fread))
  setkey(FeaMat, ID)
  Test <- FeaMat[Test_Set[, ID]]
  RFC <- readRDS(paste0("./analysis/43.ModelTraining/bandwidth/Model/NR_RF", d*1000, ".Rds"))
  data.table(bandwidth = d, Obse = Test_Set[, AA], Pred = predict(RFC, Test), Prob = apply(predict(RFC, Test, type = 'prob'), 1, max))
})

ResPer_NR <- lapply(Res, function(x) {
  x[, Obse := factor(Obse, levels = c('N', 'R'))]
  x[, Pred := factor(Pred, levels = c('N', 'R'))]
  y <- x[, confusionMatrix(Obse, Pred)]
  data.table(BandWidth = x[, unique(bandwidth)], Variable = names(c(y$overall, y$byClass)), Value = c(y$overall, y$byClass))
})
ResPer_NR <- do.call(rbind, ResPer_NR)
ResPer_NR[Variable == 'F1']


ggplot(ResPer_NR[Variable == 'F1'], aes(x = factor(BandWidth), y = Value)) + 
  geom_point()







