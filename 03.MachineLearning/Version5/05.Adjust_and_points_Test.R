setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(ggplot2)
library(parallel)
library(caret)

load(file = "./analysis/03.MachineLearning/01.data/AAInfo_RawSig.RData")

AAInfo <- AAInfo[DwellTime > 0.75 & SignalSD < 3.5 & !amino_acid %in% c("Cys", "His")]
AAInfo <- AAInfo[, .SD[sample(.N, 300), ], amino_acid]

library(doParallel)
cl <- makePSOCKcluster(20)
registerDoParallel(cl)

Adjust_Test <- lapply(c(0.1, 0.5, 1, 2, 3), function(a) {
  BinExp <- mclapply(AAInfo$ID, FUN = function(id) {
    round(density(RawSig[ID == id, pA], from = 0, to = 1, n = 100, adjust = a)$y, 3)
  }, mc.cores = 20)
  BinExp <- do.call(rbind, BinExp)
  row.names(BinExp) <- AAInfo$ID
  colnames(BinExp) <- paste0("X", sprintf("%03d", seq_len(ncol(BinExp))))
  
  Mat <- merge(as.data.table(BinExp, keep.rownames = "ID"), AAInfo[, .(ID, AreaRatio_L1, AreaRatio_L2, Outer, DwellTime, SignalSD,  Blockade, amino_acid)], by = "ID")
  Mat[, amino_acid := as.factor(amino_acid)]
  Mat <- na.omit(Mat)
  
  set.seed(3456)
  Train <- Mat[, .SD[sample(.N, 300), ], amino_acid]
  Train <- data.frame(Train[, -c(1:2)], Class = Train[[1]])
  
  fitControl <- trainControl(method = "repeatedcv",
                             number = 10,
                             repeats = 10)
  set.seed(825)
  Fit1 <- train(Class ~ ., data = Train, 
                # preProc = c("center", "scale", "YeoJohnson", "nzv"), 
                method = "rf", 
                trControl = fitControl,
                verbose = FALSE,
                ## to evaluate:
                tuneGrid = expand.grid(mtry = 36),
                # tuneLength = 10,
                metric = "Accuracy", 
                allowParallel = TRUE)
  return(Fit1)
})

stopCluster(cl)

mapply(function(x) x$results[2], Adjust_Test)
# $Accuracy
# [1] 0.7588824
# [1] 0.7668039
# [1] 0.7636471
# [1] 0.7538824
# [1] 0.7457255

Point_Test <- lapply(c(50, 100, 200, 300, 400, 500), function(p) {
  BinExp <- mclapply(AAInfo$ID, FUN = function(id) {
    round(density(RawSig[ID == id, pA], from = 0, to = 1, n = p, adjust = 0.5)$y, 3)
  }, mc.cores = 20)
  BinExp <- do.call(rbind, BinExp)
  row.names(BinExp) <- AAInfo$ID
  colnames(BinExp) <- paste0("X", sprintf("%03d", seq_len(ncol(BinExp))))
  
  Mat <- merge(as.data.table(BinExp, keep.rownames = "ID"), AAInfo[, .(ID, AreaRatio_L1, AreaRatio_L2, Outer, DwellTime, SignalSD,  Blockade, amino_acid)], by = "ID")
  Mat[, amino_acid := as.factor(amino_acid)]
  Mat <- na.omit(Mat)
  
  set.seed(3456)
  Train <- Mat[, .SD[sample(.N, 300), ], amino_acid]
  Train <- data.frame(Train[, -c(1:2)], Class = Train[[1]])
  
  fitControl <- trainControl(method = "repeatedcv",
                             number = 10,
                             repeats = 10)
  set.seed(825)
  Fit1 <- train(Class ~ ., data = Train, 
                # preProc = c("center", "scale", "YeoJohnson", "nzv"), 
                method = "rf", 
                trControl = fitControl,
                verbose = FALSE,
                ## to evaluate:
                tuneGrid = expand.grid(mtry = round(36/100*p)),
                # tuneLength = 10,
                metric = "Accuracy", 
                allowParallel = TRUE)
  return(Fit1)
})

mapply(function(x) x$results[2], Point_Test)
# $Accuracy
# [1] 0.7558627
# [1] 0.7668039
# [1] 0.774
# [1] 0.7732353
# [1] 0.7725294
# [1] 0.7721961

# $Accuracy
# [1] 0.7644706
# [1] 0.7668039
# [1] 0.7702941
# [1] 0.771098
# [1] 0.7696667
# [1] 0.7701373

# Adjust = 0.5, n = 200 似乎是最好的组合