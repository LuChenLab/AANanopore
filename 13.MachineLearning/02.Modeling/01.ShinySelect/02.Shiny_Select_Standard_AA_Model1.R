setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(caret)
Mat_List <- readRDS("./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model1/FeatueMatrix.Rds")
L1_Mat <- Mat_List[L2Ratio < .05]
setnames(L1_Mat, "A", "amino_acid")
L1_Mat$file_name <- mapply(function(x) x[1], strsplit(L1_Mat$ID, "_"))

Vali_Set <- L1_Mat[, .N, .(amino_acid, file_name)][N > 20, .SD[which.min(N)], amino_acid]
Vali_Set <- merge(Vali_Set[, .(amino_acid, file_name)], L1_Mat[, .(amino_acid, ID, file_name)], by = c("amino_acid", "file_name"))

Valid <- L1_Mat[ID %in% Vali_Set$ID]
Mat <- L1_Mat[!ID %in% Vali_Set$ID]

set.seed(3456)
Train <- Mat[, .SD[sample(.N, .N * 0.8), ], amino_acid]
set.seed(3456)
Train <- rbind(Train[amino_acid %in% Train[, .N, amino_acid][N >= 1000, amino_acid], .SD[sample(.N, 1000), ], amino_acid], 
               Train[amino_acid %in% Train[, .N, amino_acid][N < 1000, amino_acid], .SD[sample(.N, 1000, replace = T), ], amino_acid])
Test <- Mat[!ID %in% Train$ID]

Train[, ID := NULL]
Train <- data.frame(Train[, 2:609], Class = factor(Train[[1]], levels = sort(Mat[, unique(amino_acid)])))

saveRDS(Train, "./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model1/L1_Train.Rds")
saveRDS(Test, "./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model1/L1_Test.Rds")
saveRDS(Valid, "./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model1/L1_Valid.Rds")


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


Train <- readRDS("./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model1/L1_Train.Rds")
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
saveRDS(RFC, file = "./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model1/RF/RF_Fit_L1.Rds")
rm(list = c("Train", "RFC"))


Train <- readRDS("./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model1/L1_Train.Rds")
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
saveRDS(RFC, file = "./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model1/RF/RF_Fit_L1_P.Rds")
rm(list = c("Train", "RFC"))


Train <- readRDS("./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model1/L1_Train.Rds")
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
saveRDS(RFC, file = "./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model1/RF/RF_Fit_L1_X.Rds")
rm(list = c("Train", "RFC"))

