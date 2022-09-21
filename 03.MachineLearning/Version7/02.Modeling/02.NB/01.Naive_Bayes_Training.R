setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(parallel)
library(ggplot2)
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
  naive_bayes <- train(Class ~ ., 
                       data = Train, 
                       method = "naive_bayes", 
                       # tuneGrid = expand.grid(data.frame(size = 12, decay = 0.1)),
                       # tuneLength = 6,
                       trControl = fitControl)
  saveRDS(naive_bayes, file = paste0("./analysis/03.MachineLearning/Version7/02.Modeling/02.NB/01.Models/NB_Fit", i, ".Rds"))
})
stopCluster(cl)
dir.create("./analysis/03.MachineLearning/Version7/02.Modeling/02.NB/01.Models", recursive = T)
