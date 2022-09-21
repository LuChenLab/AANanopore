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

Vali_Set <- AAInfo[, .N, by = c("amino_acid", "file_name")][N > 20, .SD[which.min(N)], amino_acid]
Vali_Set <- merge(Vali_Set[, .(amino_acid, file_name)], AAInfo, by = c("amino_acid", "file_name"))

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
  round(density(RawSig[ID == id, pA], from = 0, to = 1, n = 100, adjust = 2)$y, 3)
}, mc.cores = 20)
table(mapply(length, BinExp))
BinExp <- do.call(rbind, BinExp)
row.names(BinExp) <- AAInfo1$ID
colnames(BinExp) <- paste0("X", sprintf("%03d", seq_len(ncol(BinExp))))

setkey(AAInfo1, ID)

Mat <- merge(as.data.table(BinExp, keep.rownames = "ID"), AAInfo1[row.names(BinExp), .(ID, AreaRatio_L1, AreaRatio_L2, Outer, DwellTime, SignalSD,  Blockade, amino_acid)], by = "ID")
Mat[, amino_acid := as.factor(amino_acid)]
Mat <- na.omit(Mat)

# Machine learning
# Pre-Processing of data

library(caret)

Valid <- Mat[ID %in% Vali_Set$ID]
Mat <- Mat[!ID %in% Vali_Set$ID]

set.seed(3456)
Train <- Mat[, .SD[sample(.N, 100), ], amino_acid]
Test <- Mat[!ID %in% Train$ID]

Train <- data.frame(Train[, -c(1:2)], Class = Train[[1]])

set.seed(9560)
up_Test <- upSample(x = Test[, !colnames(Test) %in% c("ID", "amino_acid"), with = F],
                    y = Test$amino_acid)                         
table(up_Test$Class)
table(Train$Class)

set.seed(9560)
up_Valid <- upSample(x = Valid[, !colnames(Valid) %in% c("ID", "amino_acid"), with = F],
                     y = Valid$amino_acid)                         
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
Fit0 <- train(Class ~ ., data = Train, 
              # preProc = c("center", "scale", "YeoJohnson", "nzv"), 
              method = "rf", 
              trControl = fitControl,
              verbose = FALSE,
              ## to evaluate:
              tuneGrid = expand.grid(mtry = 40),
              # tuneLength = 10,
              metric = "Accuracy", 
              allowParallel = TRUE)
stopCluster(cl)

varImp(Fit0)

ROC_Test <- predict(Fit0, up_Test, type = "prob")
ppm <- melt(as.data.table(ROC_Test, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
ROC_Test <- cbind(as.data.table(ROC_Test), true = up_Test$Class, ppm[, 2:3])
ROC_Test[, confusionMatrix(pred, true)]
ROC_Test[, hist(Prob, breaks = 100)]
abline(v = 0.4, col = 2)



ROC_Valid <- predict(Fit0, Valid, type = "prob")
ppm <- melt(as.data.table(ROC_Valid, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
ROC_Valid <- cbind(as.data.table(ROC_Valid), true = Valid$amino_acid, ppm[, 2:3])
ROC_Valid[, confusionMatrix(pred, true)]
ROC_Valid[Prob > 0.4, confusionMatrix(pred, true)]
ROC_Valid[, hist(Prob, breaks = 100)]
abline(v = 0.4, col = 2)
ROC_Valid[, mean(Prob > 0.4)]


res <- rbind(data.table(Dataset = "Test", ROC_Test), data.table(Dataset = "Validation", ROC_Valid))

ggplot(res, aes(Prob, colour = Dataset)) + 
  geom_line(stat = "Density", adjust = 0.5) + 
  theme_bw(base_size = 16)

ggplot(res, aes(Prob, colour = Dataset)) +
  geom_density() + 
  theme_bw(base_size = 16)

res[pred == true, Result := "Correct"]
res[pred != true, Result := "Wrong"]
ggplot(res, aes(Prob, fill = Result)) + 
  geom_histogram() + 
  theme_bw(base_size = 16) + 
  facet_wrap(~ Dataset, nrow = 2, scales = "free_y")

ggplot(res, aes(Prob, fill = Result)) + 
  geom_histogram() + 
  theme_bw(base_size = 16) + 
  facet_wrap(~ Dataset, nrow = 1, scales = "free_y") +
  theme(legend.position = "top") -> p1
p1


cutoffList <- lapply(1:50/50, function(b) {
  merge(res[, .(Recovery = mean(Prob >= b)), by = Dataset], 
        res[Prob >= b, .(Accuracy = mean(Result == "Correct")), by = Dataset])
})
cutoffList <- data.table(Cutoff = rep(1:50/50, mapply(nrow, cutoffList)), do.call(rbind, cutoffList))
cutoffList[, L := paste0(Cutoff, "(", round(Recovery * 100, 1), ", ", round(Accuracy * 100, 1), ")")]

library(ggrepel)
ggplot(cutoffList, aes(x = Recovery * 100, y = Accuracy * 100, colour = Dataset)) + 
  geom_line() + 
  theme_bw(base_size = 16) + 
  scale_x_reverse() +
  labs(x = "Recovery (%)", y = "Accuracy (%)") + 
  theme(legend.position = "top") + 
  geom_point(data = cutoffList[Cutoff %in% c(3:10/10)]) +
  geom_text_repel(data = cutoffList[Cutoff %in% c(3:10/10)], aes(label = L)) + 
  scale_color_discrete(direction = -1) -> p2

library(cowplot)
plot_grid(p1, p2, nrow = 1, rel_widths = c(2, 1.2))









