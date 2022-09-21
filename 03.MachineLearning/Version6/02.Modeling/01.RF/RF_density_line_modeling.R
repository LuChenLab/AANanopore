setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(ggplot2)
library(parallel)
library(caret)
library(ecp)

load(file = "./analysis/03.MachineLearning/01.data/Version6/AAInfo_RawSig.RData")
AAInfo <- AAInfo[!amino_acid %in% c("Cys", "Pro")]
AAInfo <- AAInfo[(DwellTime > 0.75 & SignalSD < 2.5) | amino_acid == "His"]

BinExp <- mclapply(AAInfo$ID, FUN = function(id) {
  round(density(RawSig[ID == id, pA], from = 0, to = 1, n = 200, adjust = 0.5)$y, 3)
}, mc.cores = 20)
BinExp <- do.call(rbind, BinExp)
row.names(BinExp) <- AAInfo$ID
colnames(BinExp) <- paste0("X", sprintf("%03d", seq_len(ncol(BinExp))))

Mat <- merge(as.data.table(BinExp, keep.rownames = "ID"), AAInfo[, .(ID, AreaRatio_L1, AreaRatio_L2, Outer, DwellTime, SignalSD,  Blockade, amino_acid)], by = "ID")
Mat[, amino_acid := as.factor(amino_acid)]
Mat <- na.omit(Mat)

Vali_Set <- AAInfo[, .N, by = c("amino_acid", "file_name")][N > 20, .SD[which.min(N)], amino_acid]
Vali_Set <- merge(Vali_Set[, .(amino_acid, file_name)], AAInfo, by = c("amino_acid", "file_name"))

Valid <- Mat[ID %in% Vali_Set$ID]
Mat <- Mat[!ID %in% Vali_Set$ID]

round(Mat[, .N, amino_acid][, min(N)]*0.8)

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

table(Train$Class)
table(up_Test$Class)
table(Valid$amino_acid)
table(up_Valid$Class)

save(Train, Test, up_Test, Valid, up_Valid, file = "./analysis/03.MachineLearning/01.data/Version6/Modeling_Data.RData")
save(Train, file = "./analysis/03.MachineLearning/01.data/Version6/Train_Data.RData")

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
              # tuneGrid = expand.grid(mtry = 12),
              tuneLength = 50,
              metric = "Accuracy", 
              allowParallel = TRUE)
save(Fit1, file = "./analysis/03.MachineLearning/Version6/02.Modeling/01.RF/01.Models/RF_Fit1.RData")
stopCluster(cl)

ROC_Test <- predict(Fit1, up_Test, type = "prob")
ppm <- melt(as.data.table(ROC_Test, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
ROC_Test <- cbind(as.data.table(ROC_Test), true = up_Test$Class, ppm[, 2:3])
ROC_Test[, confusionMatrix(pred, true)]
ROC_Test[, mean(pred == true)]
ROC_Test[, mean(Prob > 0.4)]
ROC_Test[Prob > 0.4, mean(pred == true)]
ROC_Test[, hist(Prob, breaks = 50)]

ROC_Valid <- predict(Fit1, Valid, type = "prob")
ppm <- melt(as.data.table(ROC_Valid, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
ROC_Valid <- cbind(as.data.table(ROC_Valid), true = Valid$amino_acid, ppm[, 2:3])
ROC_Valid[, confusionMatrix(pred, true)]
ROC_Valid[, mean(pred == true)]
ROC_Valid[, mean(Prob > 0.4)]
ROC_Valid[Prob > 0.4, mean(pred == true)]
ROC_Valid[, hist(Prob, breaks = 50)]

ROC_Valid_up <- predict(Fit1, up_Valid, type = "prob")
ppm <- melt(as.data.table(ROC_Valid_up, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
ROC_Valid_up <- cbind(as.data.table(ROC_Valid_up), true = up_Valid$Class, ppm[, 2:3])
ROC_Valid_up[, confusionMatrix(pred, true)]
ROC_Valid_up[, mean(pred == true)]
ROC_Valid_up[, mean(Prob > 0.4)]
ROC_Valid_up[Prob > 0.4, mean(pred == true)]
ROC_Valid_up[, hist(Prob, breaks = 50)]

varImp(Fit1)


par(mfrow = c(1, 2))
ROC_Test[, hist(Prob, breaks = 50)]
ROC_Valid[, hist(Prob, breaks = 50)]
dev.off()

res <- rbind(data.table(Dataset = "Test", ROC_Test), data.table(Dataset = "Validation", ROC_Valid))

ggplot(res, aes(Prob, colour = Dataset)) + 
  geom_line(stat = "Density", adjust = 0.5) + 
  theme_bw(base_size = 16)

ggplot(res, aes(Prob, fill = Dataset, colour = Dataset)) +
  geom_density(alpha = 0.4) + 
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
cutoffList <- data.table(Cutoff = rep(1:50/50, each = 2), do.call(rbind, cutoffList))
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








