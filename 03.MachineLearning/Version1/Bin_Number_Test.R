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

RawSig1 <- fread("./analysis/01.AASignalRecognition/Version9/06.MachineLearning/01.SignalPicture2/RawSignal.txt")
AAInfo <- fread("./analysis/01.AASignalRecognition/Version9/06.MachineLearning/01.SignalPicture2/AASignal.txt")


# Strategy1

answer1 <- na.omit(answer)[class != 0]
answer1[, .N, amino_acid]

AAInfo1 <- AAInfo[ID %in% answer1[, ID]]

RawSig1 <- RawSig1[ID %in% answer1[, ID]]
RawSig1 <- RawSig1[L == "U"]

RawSig1[, .N, ID][order(N)]

# x1 = matrix(RawSig1[ID == "21307011_00709", pA])
# y1 = e.divisive(X = x1, R = 1, k = 19, min.size = 3, alpha = 1)
# 
# plotSig(RawSig1[ID == id]) + 
#   geom_vline(xintercept = RawSig1[ID == id, Time][y1$estimates])

# Machine learning
# Pre-Processing of data

library(caret)
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10)

library(doParallel)
cl <- makePSOCKcluster(10)
registerDoParallel(cl)

BinTest <- list()

for(b in 9:18) {
  print(b)
  BinExp <- mclapply(AAInfo1$ID, function(id) {
    y1 = e.divisive(X = matrix(RawSig1[ID == id, pA]), sig.lvl = 0.1, R = 1, k = b, min.size = 3, alpha = 1)
    BinExp <- mapply(split(RawSig1[ID == id, pA], y1$cluster), FUN = mean2)
    return(BinExp)
  }, mc.cores = 10)
  table(mapply(length, BinExp))
  BinExp <- do.call(rbind, BinExp)
  row.names(BinExp) <- AAInfo1$ID
  colnames(BinExp) <- paste0("X", sprintf("%02d", seq_len(ncol(BinExp))))
  
  setkey(AAInfo1, ID)
  
  Mat <- merge(as.data.table(BinExp, keep.rownames = "ID"), AAInfo1[row.names(BinExp), .(ID, AreaRatio_L1, AreaRatio_L2, Outer, DwellTime, SignalMean, SignalSD,  Blockade, amino_acid)], by = "ID")
  Mat[, amino_acid := as.factor(amino_acid)]
  Mat <- na.omit(Mat)
  
  set.seed(3456)
  Train <- Mat[, .SD[sample(.N, 150), ], amino_acid]
  Test <- Mat[!ID %in% Train$ID]
  
  Train <- data.frame(Train[, -c(1:2)], Class = Train[[1]])
  
  set.seed(9560)
  up_Test <- upSample(x = Test[, !colnames(Test) %in% c("ID", "amino_acid"), with = F],
                      y = Test$amino_acid)                         
  table(up_Test$Class)
  table(Train$Class)
  
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
  
  ROC_Tab <- predict(Fit1, up_Test, type = "prob")
  ppm <- melt(as.data.table(ROC_Tab, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
  ROC_Tab <- cbind(as.data.table(ROC_Tab), true = up_Test$Class, ppm[, 2:3])
  cM <- ROC_Tab[, confusionMatrix(pred, true)]
  
  BinTest[[b]] <- data.table(Bin = b, X = paste0(Fit1$coefnames, collapse = ","), 
                             TrainAccuracy = max(Fit1$results[, 2]), 
                             TestAccuracy = as.numeric(cM$overall[1]))
}

stopCluster(cl)

do.call(rbind, BinTest)
# Bin                                                                                                                                                  X TrainAccuracy TestAccuracy
# 1:   9                                     X01,X02,X03,X04,X05,X06,X07,X08,X09,X10,AreaRatio_L1,AreaRatio_L2,Outer,DwellTime,SignalMean,SignalSD,Blockade     0.8997895    0.8969802
# 2:  10                                 X01,X02,X03,X04,X05,X06,X07,X08,X09,X10,X11,AreaRatio_L1,AreaRatio_L2,Outer,DwellTime,SignalMean,SignalSD,Blockade     0.8461053    0.8421053
# 3:  11                             X01,X02,X03,X04,X05,X06,X07,X08,X09,X10,X11,X12,AreaRatio_L1,AreaRatio_L2,Outer,DwellTime,SignalMean,SignalSD,Blockade     0.8452632    0.8481450
# 4:  12                         X01,X02,X03,X04,X05,X06,X07,X08,X09,X10,X11,X12,X13,AreaRatio_L1,AreaRatio_L2,Outer,DwellTime,SignalMean,SignalSD,Blockade     0.8416842    0.8500431
# 5:  13                     X01,X02,X03,X04,X05,X06,X07,X08,X09,X10,X11,X12,X13,X14,AreaRatio_L1,AreaRatio_L2,Outer,DwellTime,SignalMean,SignalSD,Blockade     0.8435789    0.8479724
# 6:  14                 X01,X02,X03,X04,X05,X06,X07,X08,X09,X10,X11,X12,X13,X14,X15,AreaRatio_L1,AreaRatio_L2,Outer,DwellTime,SignalMean,SignalSD,Blockade     0.8451579    0.8460742
# 7:  15             X01,X02,X03,X04,X05,X06,X07,X08,X09,X10,X11,X12,X13,X14,X15,X16,AreaRatio_L1,AreaRatio_L2,Outer,DwellTime,SignalMean,SignalSD,Blockade     0.8430175    0.8476273
# 8:  16         X01,X02,X03,X04,X05,X06,X07,X08,X09,X10,X11,X12,X13,X14,X15,X16,X17,AreaRatio_L1,AreaRatio_L2,Outer,DwellTime,SignalMean,SignalSD,Blockade     0.8441404    0.8431406
# 9:  17     X01,X02,X03,X04,X05,X06,X07,X08,X09,X10,X11,X12,X13,X14,X15,X16,X17,X18,AreaRatio_L1,AreaRatio_L2,Outer,DwellTime,SignalMean,SignalSD,Blockade     0.8410526    0.8479724
# 10:  18 X01,X02,X03,X04,X05,X06,X07,X08,X09,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,AreaRatio_L1,AreaRatio_L2,Outer,DwellTime,SignalMean,SignalSD,Blockade     0.8416842    0.8455565


