setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(ggplot2)
library(parallel)
library(caret)

c("Blockade", "AllTime", "DwellTime", "SignalSD")

fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10)

library(doParallel)
cl <- makePSOCKcluster(20)
registerDoParallel(cl)

lapply(1:6, function(i) {
  print(i)
  load(file = paste0("./analysis/03.MachineLearning/01.data/Version7/Train_Data_", i, ".RData"))
  set.seed(825)
  SVM <- train(Class ~ ., 
               data = Train, 
               method = "lssvmRadial", 
               trControl = fitControl,
               verbose = FALSE)
  saveRDS(SVM, file = paste0("./analysis/03.MachineLearning/Version7/02.Modeling/03.SVM/01.Models/SVM_Fit", i, ".Rds"))
})
stopCluster(cl)
















