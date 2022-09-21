setwd("/mnt/raid61/Personal_data/tangchao/AANanopore/")
library(data.table)
library(ecp)
library(ggplot2)
library(parallel)

mean2 <- function(x, q1 = 0.25, q2 = 0.75) {
  q <- quantile(x, c(q1, q2))
  mean(x[x >= min(q) & x <= max(q)])
}


plotSig <- function(x) {
  ggplot() + 
    geom_line(data = x, mapping = aes(x = Time, y = pA)) + 
    geom_line(data = x[L == "U"], mapping = aes(x = Time, y = pA), colour = "red", size = 1.1) + 
    theme_minimal(base_size = 15)
}

RawSig <- fread("./analysis/01.AASignalRecognition/Version9/06.MachineLearning/01.SignalPicture2/RawSignal.txt")
AAInfo <- fread("./analysis/01.AASignalRecognition/Version9/06.MachineLearning/01.SignalPicture2/AASignal.txt")
AAInfo <- na.omit(AAInfo)

RawSig <- RawSig[L == "U"]

# Strategy1

RawSig[, .N, ID][order(N)]

# x1 = matrix(RawSig[ID == "21307011_00709", pA])
# y1 = e.divisive(X = x1, R = 1, k = 19, min.size = 3, alpha = 1)
# 
# plotSig(RawSig[ID == id]) + 
#   geom_vline(xintercept = RawSig[ID == id, Time][y1$estimates])


BinExp <- mclapply(AAInfo$ID, function(id) {
  y1 <- tryCatch(e.divisive(X = matrix(RawSig[ID == id, pA]), sig.lvl = 0.1, R = 1, k = 15, min.size = 4, alpha = 1), error = function(e) NA)
  BinExp <- tryCatch(mapply(split(RawSig[ID == id, pA], y1$cluster), FUN = mean2), error = function(e) NA)
  return(BinExp)
}, mc.cores = 20)
names(BinExp) <- AAInfo$ID
BinExp <- BinExp[mapply(length, BinExp) == 16]
table(mapply(length, BinExp))

BinExp <- do.call(rbind, BinExp)
colnames(BinExp) <- paste0("X", sprintf("%02d", seq_len(ncol(BinExp))))

Mat <- merge(as.data.table(BinExp, keep.rownames = "ID"), AAInfo[, .(ID, AreaRatio_L1, AreaRatio_L2, Outer, DwellTime, SignalMean, SignalSD,  Blockade, amino_acid)], by = "ID")
Mat[, amino_acid := as.factor(amino_acid)]
Mat <- na.omit(Mat)

# Machine learning
# Pre-Processing of data

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
              tuneGrid = expand.grid(mtry = 5:12),
              # tuneLength = 10,
              metric = "Accuracy", 
              allowParallel = TRUE)

set.seed(825)
Fit2 <- train(Class ~ ., data = Train, 
              preProc = c("center", "scale", "YeoJohnson", "nzv"),
              method = "rf", 
              trControl = fitControl,
              verbose = FALSE,
              ## to evaluate:
              tuneGrid = expand.grid(mtry = 5:18),
              # tuneLength = 10,
              metric = "Accuracy", 
              allowParallel = TRUE)
stopCluster(cl)


ROC_Tab1 <- predict(Fit1, up_Test, type = "prob")
ppm <- melt(as.data.table(ROC_Tab1, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
ROC_Tab1 <- cbind(as.data.table(ROC_Tab1), true = up_Test$Class, ppm[, 2:3])
ROC_Tab1[, confusionMatrix(pred, true)]


ROC_Tab2 <- predict(Fit2, up_Test, type = "prob")
ppm <- melt(as.data.table(ROC_Tab2, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
ROC_Tab2 <- cbind(as.data.table(ROC_Tab2), true = up_Test$Class, ppm[, 2:3])
ROC_Tab2[, confusionMatrix(pred, true)]



