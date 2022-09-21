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
  pcaNNet <- train(Class ~ ., 
                   data = Train, 
                   method = "pcaNNet", 
                   # tuneGrid = expand.grid(data.frame(size = 12, decay = 0.1)),
                   tuneLength = 6,
                   trControl = fitControl)
  saveRDS(pcaNNet, file = paste0("./analysis/03.MachineLearning/Version7/02.Modeling/04.pcaNNet/01.Models/pcaNNet_Fit", i, ".Rds"))
})
stopCluster(cl)
dir.create("./analysis/03.MachineLearning/Version7/02.Modeling/04.pcaNNet/01.Models", recursive = T)











