setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(parallel)

files2 <- list.files("./analysis/11.SignalIdentification/Jan07/AASignal", pattern = "_Selected.txt", full.names = TRUE)
List <- lapply(files2, fread)
names(List) <- gsub("RawSignal_", "", gsub("_Selected.txt", "", basename(files2)))

Mat_List <- mclapply(seq_along(List), function(i) {
  print(i)
  Current <- readRDS(paste0("./analysis/11.SignalIdentification/Dec27/SignalCurrent_", names(List)[i], ".Rds"))
  Mat <- List[[i]]
  Feature_Mat <- lapply(Mat$ID, function(j) {
    D <- density(Current[ID == j, Current], from = 0, to = 1, n = 500, adjust = 0.5)$y
    pA <- sort(Current[ID == j, Current])
    pA <- pA[round(seq(1, length(pA), length.out = 100))]
    D <- c(D, pA)
    cbind(t(data.frame(round(D, 4), row.names = c(paste0("X", sprintf("%03d", 1:500)), paste0("P", sprintf("%03d", 1:100))))), 
          Mat[ID == j, .(DeltaMean, StageSD, CurrentSD, Segments, SignalCurrentPercent, L2Ratio, DwellTime, Blockade, A, ID)])
  })
  do.call(rbind, Feature_Mat)
}, mc.cores = 40)
Mat_List <- do.call(rbind, Mat_List)

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

saveRDS(Train, "./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model0/L1_Train.Rds")
saveRDS(Test, "./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model0/L1_Test.Rds")
saveRDS(Valid, "./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model0/L1_Valid.Rds")


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


Train <- readRDS("./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model0/L1_Train.Rds")
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
saveRDS(RFC, file = "./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model0/RF/RF_Fit_L1.Rds")
rm(list = c("Train", "RFC"))


Train <- readRDS("./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model0/L1_Train.Rds")
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
saveRDS(RFC, file = "./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model0/RF/RF_Fit_L1_P.Rds")
rm(list = c("Train", "RFC"))


Train <- readRDS("./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model0/L1_Train.Rds")
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
saveRDS(RFC, file = "./analysis/13.MachineLearning/01.DataPreparation/01.ShinySelect/Model0/RF/RF_Fit_L1_X.Rds")
rm(list = c("Train", "RFC"))






