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


answer <- as.data.table(openxlsx::read.xlsx("./analysis/03.MachineLearning/01.data/AASignal_20211206.xlsx"))
colnames(answer) <- c("amino_acid", "ID", "class")
answer[class != 0, .N, amino_acid]

RawSig <- fread("./analysis/01.AASignalRecognition/Version9/06.MachineLearning/01.SignalPicture2/RawSignal.txt")
AAInfo <- fread("./analysis/01.AASignalRecognition/Version9/06.MachineLearning/01.SignalPicture2/AASignal.txt")


# Strategy1

answer1 <- na.omit(answer)[class != 0]
answer1[, .N, amino_acid]

AAInfo1 <- AAInfo[ID %in% answer1[, ID]]

RawSig1 <- RawSig[ID %in% answer1[, ID]]
RawSig1 <- RawSig1[L == "U"]

RawSig1[, .N, ID][order(N)]

# x1 = matrix(RawSig1[ID == "21307011_00709", pA])
# y1 = e.divisive(X = x1, R = 1, k = 19, min.size = 3, alpha = 1)
# 
# plotSig(RawSig1[ID == id]) + 
#   geom_vline(xintercept = RawSig1[ID == id, Time][y1$estimates])


BinExp <- mclapply(AAInfo1$ID, function(id) {
  y1 = e.divisive(X = matrix(RawSig1[ID == id, pA]), sig.lvl = 0.1, R = 1, k = 9, min.size = 3, alpha = 1)
  BinExp <- mapply(split(RawSig1[ID == id, pA], y1$cluster), FUN = mean2)
  return(BinExp)
}, mc.cores = 20)
table(mapply(length, BinExp))
BinExp <- do.call(rbind, BinExp)
row.names(BinExp) <- AAInfo1$ID
colnames(BinExp) <- paste0("X", sprintf("%02d", seq_len(ncol(BinExp))))

setkey(AAInfo1, ID)

Mat <- merge(as.data.table(BinExp, keep.rownames = "ID"), AAInfo1[row.names(BinExp), .(ID, AreaRatio_L1, AreaRatio_L2, Outer, DwellTime, SignalMean, SignalSD,  Blockade, amino_acid)], by = "ID")
Mat[, amino_acid := as.factor(amino_acid)]
Mat <- na.omit(Mat)

# Machine learning
# Pre-Processing of data

library(caret)

set.seed(3456)
Train <- Mat[, .SD[sample(.N, 150), ], amino_acid]
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
              tuneGrid = expand.grid(mtry = 5:15),
              # tuneLength = 10,
              metric = "Accuracy", 
              allowParallel = TRUE)
stopCluster(cl)


ROC_Tab <- predict(Fit1, up_Test, type = "prob")
ppm <- melt(as.data.table(ROC_Tab, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
ROC_Tab <- cbind(as.data.table(ROC_Tab), true = up_Test$Class, ppm[, 2:3])

ROC_Tab[, confusionMatrix(pred, true)]

