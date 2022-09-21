setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(ggplot2)
library(parallel)
library(caret)

c("Blockade", "AllTime", "DwellTime", "SignalSD")

load(file = "./analysis/03.MachineLearning/01.data/Version6/AAInfo_RawSig.RData")
AAInfo <- AAInfo[!amino_acid %in% c("Cys", "Pro")]
AAInfo$AllTime <- mapply(AAInfo$ID, FUN = function(x) {
  RawSig[ID == x, diff(range(Time))] * 1000
})

fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10)

library(doParallel)
cl <- makePSOCKcluster(20)
registerDoParallel(cl)

lapply(1:6, function(i) {
  AAInfoi <- AAInfo[AllTime > i & SignalSD < 4]
  
  BinExp <- mclapply(AAInfoi$ID, FUN = function(id) {
    round(density(RawSig[ID == id, pA], from = 0, to = 1, n = 200, adjust = 0.5)$y, 3)
  }, mc.cores = 20)
  BinExp <- do.call(rbind, BinExp)
  row.names(BinExp) <- AAInfoi$ID
  colnames(BinExp) <- paste0("X", sprintf("%03d", seq_len(ncol(BinExp))))
  
  Mat <- merge(as.data.table(BinExp, keep.rownames = "ID"), AAInfoi[, .(ID, AllTime, DwellTime, SignalSD, Blockade, amino_acid)], by = "ID")
  Mat[, amino_acid := as.factor(amino_acid)]
  Mat <- na.omit(Mat)
  
  Vali_Set <- AAInfoi[, .N, by = c("amino_acid", "file_name")][N > 20, .SD[which.min(N)], amino_acid]
  Vali_Set <- merge(Vali_Set[, .(amino_acid, file_name)], AAInfoi, by = c("amino_acid", "file_name"))
  
  Valid <- Mat[ID %in% Vali_Set$ID]
  Mat <- Mat[!ID %in% Vali_Set$ID]
  
  set.seed(3456)
  Train <- Mat[, .SD[sample(.N, round(Mat[, .N, amino_acid][, min(N)]*0.8)), ], amino_acid]
  Test <- Mat[!ID %in% Train$ID]
  
  Train <- data.frame(Train[, -c(1:2)], Class = Train[[1]])
  
  set.seed(9560)
  up_Test <- upSample(x = Test[, !colnames(Test) %in% c("ID", "amino_acid"), with = F],
                      y = Test$amino_acid)                         
  
  set.seed(9560)
  up_Valid <- upSample(x = Valid[, !colnames(Valid) %in% c("ID", "amino_acid"), with = F],
                       y = Valid$amino_acid)                         
  
  save(Train, Test, up_Test, Valid, up_Valid, file = paste0("./analysis/03.MachineLearning/01.data/Version7/Modeling_Data_", i, ".RData"))
  save(Train, file = paste0("./analysis/03.MachineLearning/01.data/Version7/Train_Data_", i, ".RData"))
  
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
  saveRDS(RFC, file = paste0("./analysis/03.MachineLearning/Version7/02.Modeling/01.RF/01.Models/RF_Fit", i, ".Rds"))
})
stopCluster(cl)
