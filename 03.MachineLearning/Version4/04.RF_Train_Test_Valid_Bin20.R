setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(ggplot2)
library(parallel)
library(ecp)

mean2 <- function(x, q1 = 0.25, q2 = 0.75) {
  q <- quantile(x, c(q1, q2))
  mean(x[x >= min(q) & x <= max(q)])
}

load(file = "./analysis/03.MachineLearning/01.data/AAInfo_RawSig.RData")

BinExp <- mclapply(AAInfo$ID, function(id) {
  y1 <- tryCatch(e.divisive(X = matrix(RawSig[ID == id, pA]), sig.lvl = 0.1, R = 1, k = 19, min.size = 4, alpha = 1), error = function(e) NA)
  BinExp <- tryCatch(mapply(split(RawSig[ID == id, pA], y1$cluster), FUN = mean2), error = function(e) NA)
  return(BinExp)
}, mc.cores = 10)
table(mapply(length, BinExp))
names(BinExp) <- AAInfo$ID
BinExp <- BinExp[mapply(length, BinExp) == max(mapply(length, BinExp))]
table(mapply(length, BinExp))

BinExp <- do.call(rbind, BinExp)
colnames(BinExp) <- paste0("X", sprintf("%02d", seq_len(ncol(BinExp))))

save(AAInfo, RawSig, BinExp, file = "./analysis/03.MachineLearning/01.data/AAInfo_RawSig_BinExp2.RData")


Mat <- merge(as.data.table(BinExp, keep.rownames = "ID"), AAInfo[, .(ID, AreaRatio_L1, AreaRatio_L2, Outer, DwellTime, SignalMean, SignalSD,  Blockade, amino_acid)], by = "ID")
Mat[, amino_acid := as.factor(amino_acid)]
Mat <- na.omit(Mat)

Vali_Set <- AAInfo[, .N, by = c("amino_acid", "file_name")][N > 20, .SD[which.min(N)], amino_acid]
Vali_Set <- merge(Vali_Set[, .(amino_acid, file_name)], AAInfo, by = c("amino_acid", "file_name"))

library(caret)

Valid <- Mat[ID %in% Vali_Set$ID]
Mat <- Mat[!ID %in% Vali_Set$ID]

round(Mat[, .N, amino_acid][, min(N)]*0.75)

set.seed(3456)
Train <- Mat[, .SD[sample(.N, round(Mat[, .N, amino_acid][, min(N)]*0.75)), ], amino_acid]
Test <- Mat[!ID %in% Train$ID]

Train <- data.frame(Train[, -c(1:2)], Class = Train[[1]])

set.seed(9560)
up_Test <- upSample(x = Test[, !colnames(Test) %in% c("ID", "amino_acid"), with = F],
                    y = Test$amino_acid)                         

set.seed(9560)
up_Valid <- upSample(x = Valid[, !colnames(Valid) %in% c("ID", "amino_acid"), with = F],
                     y = Valid$amino_acid)                         

table(Train$Class)
table(up_Test$Class)
table(Valid$amino_acid)
table(up_Valid$Class)

fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10)



library(doParallel)
cl <- makePSOCKcluster(20)
registerDoParallel(cl)

# library(doSNOW)
# cl <- makeSOCKcluster(6)
# registerDoSNOW(cl)

# rfGrid <-  expand.grid(mtry = 2)

set.seed(825)
Fit1 <- train(Class ~ ., data = Train, 
              # preProc = c("center", "scale", "YeoJohnson", "nzv"), 
              method = "rf", 
              trControl = fitControl,
              verbose = FALSE,
              ## to evaluate:
              tuneGrid = expand.grid(mtry = 12),
              # tuneLength = 10,
              metric = "Accuracy", 
              allowParallel = TRUE)

Fit2 <- train(Class ~ ., data = Train, 
              # preProc = c("center", "scale", "YeoJohnson", "nzv"), 
              method = "rf", 
              trControl = trainControl(method = "LOOCV"),
              verbose = FALSE,
              ## to evaluate:
              tuneGrid = expand.grid(mtry = 12),
              # tuneLength = 10,
              metric = "Accuracy", 
              allowParallel = TRUE)
stopCluster(cl)

ROC_Test <- predict(Fit1, up_Test, type = "prob")
ppm <- melt(as.data.table(ROC_Test, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
ROC_Test <- cbind(as.data.table(ROC_Test), true = up_Test$Class, ppm[, 2:3])
ROC_Test[, confusionMatrix(pred, true)]
ROC_Test[, mean(pred == true)]
ROC_Test[, mean(Prob > 0.4)]
ROC_Test[Prob > 0.4, mean(pred == true)]

ROC_Test2 <- predict(Fit2, up_Test, type = "prob")
ppm <- melt(as.data.table(ROC_Test2, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
ROC_Test2 <- cbind(as.data.table(ROC_Test2), true = up_Test$Class, ppm[, 2:3])
ROC_Test2[, confusionMatrix(pred, true)]
ROC_Test2[, mean(Prob > 0.4)]
ROC_Test2[Prob > 0.4, mean(pred == true)]


ROC_Valid <- predict(Fit1, Valid, type = "prob")
ppm <- melt(as.data.table(ROC_Valid, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
ROC_Valid <- cbind(as.data.table(ROC_Valid), true = Valid$amino_acid, ppm[, 2:3])
ROC_Valid[, confusionMatrix(pred, true)]
ROC_Valid[, mean(pred == true)]
ROC_Valid[, mean(Prob > 0.4)]
ROC_Valid[Prob > 0.4, mean(pred == true)]

ROC_Valid2 <- predict(Fit2, Valid, type = "prob")
ppm <- melt(as.data.table(ROC_Valid2, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
ROC_Valid2 <- cbind(as.data.table(ROC_Valid2), true = Valid$amino_acid, ppm[, 2:3])
ROC_Valid2[, confusionMatrix(pred, true)]
ROC_Valid2[, mean(Prob > 0.4)]
ROC_Valid2[Prob > 0.4, mean(pred == true)]

ROC_Valid_up <- predict(Fit1, up_Valid, type = "prob")
ppm <- melt(as.data.table(ROC_Valid_up, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
ROC_Valid_up <- cbind(as.data.table(ROC_Valid_up), true = up_Valid$Class, ppm[, 2:3])
ROC_Valid_up[, confusionMatrix(pred, true)]
ROC_Valid_up[, mean(pred == true)]
ROC_Valid_up[, mean(Prob > 0.4)]
ROC_Valid_up[Prob > 0.4, mean(pred == true)]

ROC_Valid_up2 <- predict(Fit2, up_Valid, type = "prob")
ppm <- melt(as.data.table(ROC_Valid_up2, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
ROC_Valid_up2 <- cbind(as.data.table(ROC_Valid_up2), true = up_Valid$Class, ppm[, 2:3])
ROC_Valid_up2[, confusionMatrix(pred, true)]
ROC_Valid_up2[, mean(pred == true)]
ROC_Valid_up2[, mean(Prob > 0.4)]
ROC_Valid_up2[Prob > 0.4, mean(pred == true)]












