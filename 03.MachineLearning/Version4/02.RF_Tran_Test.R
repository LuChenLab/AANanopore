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


load(file = "./analysis/03.MachineLearning/01.data/AAInfo_RawSig_BinExp.RData")

# Strategy1

BinExp <- mclapply(AAInfo$ID, function(id) {
  y1 <- tryCatch(e.divisive(X = matrix(RawSig[ID == id, pA]), sig.lvl = 0.1, R = 1, k = 10, min.size = 4, alpha = 1), error = function(e) NA)
  BinExp <- tryCatch(mapply(split(RawSig[ID == id, pA], y1$cluster), FUN = mean2), error = function(e) NA)
  return(BinExp)
}, mc.cores = 10)
names(BinExp) <- AAInfo$ID
BinExp <- BinExp[mapply(length, BinExp) == max(mapply(length, BinExp))]
table(mapply(length, BinExp))

BinExp <- do.call(rbind, BinExp)
colnames(BinExp) <- paste0("X", sprintf("%02d", seq_len(ncol(BinExp))))

save(AAInfo, RawSig, BinExp, file = "./analysis/03.MachineLearning/01.data/AAInfo_RawSig_BinExp.RData")

Mat <- merge(as.data.table(BinExp, keep.rownames = "ID"), AAInfo[, .(ID, AreaRatio_L1, AreaRatio_L2, Outer, DwellTime, SignalMean, SignalSD,  Blockade, amino_acid)], by = "ID")
Mat[, amino_acid := as.factor(amino_acid)]
Mat <- na.omit(Mat)


library(caret)

set.seed(3456)
Train <- Mat[, .SD[sample(.N, round(Mat[, .N, amino_acid][, min(N)]*0.7)), ], amino_acid]
Test <- Mat[!ID %in% Train$ID]

Train <- data.frame(Train[, -c(1:2)], Class = Train[[1]])

set.seed(9560)
up_Test <- upSample(x = Test[, !colnames(Test) %in% c("ID", "amino_acid"), with = F],
                    y = Test$amino_acid)                         
table(up_Test$Class)
table(Train$Class)


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
              tuneGrid = expand.grid(mtry = 11),
              # tuneLength = 10,
              metric = "Accuracy", 
              allowParallel = TRUE)

stopCluster(cl)


ROC_Tab1 <- predict(Fit1, up_Test, type = "prob")
ppm <- melt(as.data.table(ROC_Tab1, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
ROC_Tab1 <- cbind(as.data.table(ROC_Tab1), true = up_Test$Class, ppm[, 2:3])
ROC_Tab1[, confusionMatrix(pred, true)]
ROC_Tab1[, mean(Prob > 0.4)]
ROC_Tab1[Prob > 0.4, mean(pred == true)]


