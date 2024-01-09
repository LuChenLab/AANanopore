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


Train <- readRDS("./analysis/13.MachineLearning/01.DataPreparation/Jan10/Model0/L1_Train.Rds")
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
saveRDS(RFC, file = "./analysis/13.MachineLearning/01.DataPreparation/Jan10/Model0/RF/RF_Fit_L1.Rds")
rm(list = c("Train", "RFC"))


Train <- readRDS("./analysis/13.MachineLearning/01.DataPreparation/Jan10/Model0/L1_Train.Rds")
Train <- subset(Train, select = !colnames(Train) %in% paste0("X", sprintf("%03d", 1:500)))
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
saveRDS(RFC, file = "./analysis/13.MachineLearning/01.DataPreparation/Jan10/Model0/RF/RF_Fit_L1_P.Rds")
rm(list = c("Train", "RFC"))



Train <- readRDS("./analysis/13.MachineLearning/01.DataPreparation/Jan10/Model0/L1_Train.Rds")
Train <- subset(Train, select = !colnames(Train) %in% paste0("P", sprintf("%03d", 1:100)))
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
saveRDS(RFC, file = "./analysis/13.MachineLearning/01.DataPreparation/Jan10/Model0/RF/RF_Fit_L1_X.Rds")
rm(list = c("Train", "RFC"))

















Train <- readRDS("./analysis/13.MachineLearning/01.DataPreparation/Jan10/Model0/L2_Train.Rds")
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
saveRDS(RFC, file = "./analysis/13.MachineLearning/01.DataPreparation/Jan10/Model0/RF/RF_Fit_L2.Rds")
rm(list = c("Train", "RFC"))


Train <- readRDS("./analysis/13.MachineLearning/01.DataPreparation/Jan10/Model0/L2_Train.Rds")
Train <- subset(Train, select = !colnames(Train) %in% paste0("X", sprintf("%03d", 1:500)))
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
saveRDS(RFC, file = "./analysis/13.MachineLearning/01.DataPreparation/Jan10/Model0/RF/RF_Fit_L2_P.Rds")
rm(list = c("Train", "RFC"))



Train <- readRDS("./analysis/13.MachineLearning/01.DataPreparation/Jan10/Model0/L2_Train.Rds")
Train <- subset(Train, select = !colnames(Train) %in% paste0("P", sprintf("%03d", 1:100)))
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
saveRDS(RFC, file = "./analysis/13.MachineLearning/01.DataPreparation/Jan10/Model0/RF/RF_Fit_L2_X.Rds")
rm(list = c("Train", "RFC"))












