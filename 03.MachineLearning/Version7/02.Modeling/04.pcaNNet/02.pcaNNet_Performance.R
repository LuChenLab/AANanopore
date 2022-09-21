setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(ggplot2)
library(parallel)
library(caret)

Modeling_Datas <- list.files("./analysis/03.MachineLearning/01.data/Version7", "Modeling_Data", full.names = T)
Models <- list.files("./analysis/03.MachineLearning/Version7/02.Modeling/04.pcaNNet/01.Models", "Rds", full.names = T)

Preds <- mclapply(1:6, function(i) {
  load(Modeling_Datas[i])
  Fit <- readRDS(Models[i])
  
  ROC_Test <- predict(Fit, Test, type = "prob")
  ppm <- melt(as.data.table(ROC_Test, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
  ROC_Test <- cbind(as.data.table(ROC_Test), true = Test$amino_acid, ppm[, 2:3])
  
  ROC_Test_up <- predict(Fit, up_Test, type = "prob")
  ppm <- melt(as.data.table(ROC_Test_up, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
  ROC_Test_up <- cbind(as.data.table(ROC_Test_up), true = up_Test$Class, ppm[, 2:3])
  
  ROC_Valid <- predict(Fit, Valid, type = "prob")
  ppm <- melt(as.data.table(ROC_Valid, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
  ROC_Valid <- cbind(as.data.table(ROC_Valid), true = Valid$amino_acid, ppm[, 2:3])
  
  ROC_Valid_up <- predict(Fit, up_Valid, type = "prob")
  ppm <- melt(as.data.table(ROC_Valid_up, keep.rownames = "ID"), value.name = "Prob", variable.name = "pred")[, .SD[which.max(Prob), ], ID]
  ROC_Valid_up <- cbind(as.data.table(ROC_Valid_up), true = up_Valid$Class, ppm[, 2:3])
  
  res <- rbind(data.table(Dataset = "Test", ROC_Test), 
               data.table(Dataset = "TestUp", ROC_Test_up),
               data.table(Dataset = "Validation", ROC_Valid),
               data.table(Dataset = "ValidationUp", ROC_Valid_up))
  data.table(Model = paste0("Fit", i), res)
}, mc.cores = 6)


cutoffList <- mclapply(Preds, function(res) {
  res[pred == true, Result := "Correct"]
  res[pred != true, Result := "Wrong"]
  cutoffList <- lapply(1:50/50, function(b) {
    merge(res[, .(Recovery = mean(Prob >= b)), by = Dataset], 
          res[Prob >= b, .(Accuracy = mean(Result == "Correct")), by = Dataset])
  })
  cutoffList <- data.table(Cutoff = rep(1:50/50, mapply(nrow, cutoffList)), do.call(rbind, cutoffList))
  cutoffList[, L := paste0(Cutoff, "(", round(Recovery * 100, 1), ", ", round(Accuracy * 100, 1), ")")]
  return(cutoffList)
}, mc.cores = 6)

cutoffList <- data.table(Model = rep(paste0("pcaNNet", 1:6), mapply(nrow, cutoffList)), do.call(rbind, cutoffList))

library(ggrepel)
ggplot(cutoffList[Dataset == "Test"], aes(x = Recovery * 100, y = Accuracy * 100, colour = Model)) + 
  geom_line() + 
  theme_bw(base_size = 16) + 
  scale_x_reverse() +
  labs(x = "Recovery (%)", y = "Accuracy (%)", title = "Test set") + 
  theme(legend.position = "top") + 
  geom_point(data = cutoffList[Dataset == "Test" & Cutoff %in% c(c(3,5,7,10)/10)]) +
  geom_text_repel(data = cutoffList[Dataset == "Test" & Cutoff %in% c(c(3,5,7,10)/10)], aes(label = L)) + 
  scale_color_brewer(palette = "Dark2")


ggplot(cutoffList[Dataset == "Validation"], aes(x = Recovery * 100, y = Accuracy * 100, colour = Model)) + 
  geom_line() + 
  theme_bw(base_size = 16) + 
  scale_x_reverse() +
  labs(x = "Recovery (%)", y = "Accuracy (%)", title = "Validation set") + 
  theme(legend.position = "top") + 
  geom_point(data = cutoffList[Dataset == "Validation" & Cutoff %in% c(c(3,5,7,10)/10)]) +
  geom_text_repel(data = cutoffList[Dataset == "Validation" & Cutoff %in% c(c(3,5,7,10)/10)], aes(label = L)) + 
  scale_color_brewer(palette = "Dark2")





openxlsx::write.xlsx(cutoffList, "analysis/03.MachineLearning/Version7/02.Modeling/04.pcaNNet/02.CompareModels/cutoffList.xlsx")
saveRDS(Preds, file = "analysis/03.MachineLearning/Version7/02.Modeling/04.pcaNNet/02.CompareModels/Preds.Rds")

cutoffList <- as.data.table(read.xlsx("analysis/03.MachineLearning/Version7/02.Modeling/04.pcaNNet/02.CompareModels/cutoffList.xlsx"))
Preds <- readRDS("analysis/03.MachineLearning/Version7/02.Modeling/04.pcaNNet/02.CompareModels/Preds.Rds")

library(ggrepel)

ggplot(cutoffList[Dataset %in% c("Test", "Validation")], aes(x = Recovery * 100, y = Accuracy * 100, colour = Model)) + 
  geom_line() + 
  theme_bw(base_size = 16) + 
  scale_x_reverse() +
  labs(x = "Recovery (%)", y = "Accuracy (%)") + 
  theme(legend.position = "top") + 
  geom_point(data = cutoffList[Dataset %in% c("Test", "Validation") & Cutoff %in% c(c(7)/10)]) +
  geom_text_repel(data = cutoffList[Dataset %in% c("Test", "Validation") & Cutoff %in% c(c(7)/10)], aes(label = L)) + 
  scale_color_brewer(palette = "Dark2") + 
  facet_wrap(~ Dataset) + 
  guides(colour = guide_legend(nrow = 1))
ggsave("analysis/03.MachineLearning/Version7/02.Modeling/04.pcaNNet/02.CompareModels/Cutoff_Recovery_Accuracy.pdf", width = 9, height = 6)

ggplot(SignalsNumber[, mean(Train), Model], aes(x = Model, y = V1, fill = Model)) + 
  geom_col() + 
  geom_text(aes(label = V1, y = V1 + 20)) + 
  scale_fill_brewer(palette = "Dark2") + 
  labs(y = "# training set") + 
  theme_bw(base_size = 15) -> p1
p1 <- p1 + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = "none")

ggplot(SignalsNumber, aes(x = Model, y = Test)) + 
  geom_boxplot(aes(colour = Model)) + 
  scale_colour_brewer(palette = "Dark2") + 
  labs(y = "# test set") + 
  theme_bw(base_size = 15) + 
  geom_text(data = SignalsNumber[, median(Test), Model], 
            mapping = aes(y = V1 + 200, label = V1)) -> p2
p2 <- p2 + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = "none")

ggplot(SignalsNumber, aes(x = Model, y = Validation)) + 
  geom_boxplot(aes(colour = Model)) + 
  scale_colour_brewer(palette = "Dark2") + 
  labs(y = "# validation set") + 
  guides(colour = "none") +
  theme_bw(base_size = 15) + 
  geom_text(data = SignalsNumber[, median(Validation), Model], 
            mapping = aes(y = V1 + 30, label = V1)) -> p3


library(patchwork)
p1 / p2 / p3
ggsave("analysis/03.MachineLearning/Version7/02.Modeling/04.pcaNNet/02.CompareModels/Signal_Number_of_Train_Test_Validation.pdf", width = 5, height = 6)



Preds <- do.call(rbind, Preds)
Preds[pred == true, Result := "Correct"]
Preds[pred != true, Result := "Wrong"]
Preds[, Model := gsub("Fit", "pcaNNet", Model)]


ggplot(Preds[Dataset %in% c("Test", "Validation")], aes(Prob, fill = Result)) +
  geom_histogram() +
  theme_bw(base_size = 16) +
  facet_wrap(~ Dataset + Model, nrow = 2, scales = "free_y") +
  theme(legend.position = "top") + 
  labs(x = "Probability", y = "Count")

ggsave("analysis/03.MachineLearning/Version7/02.Modeling/04.pcaNNet/02.CompareModels/Prediction_Probability_Distribution.pdf", width = 16, height = 6)




ggplot(Preds[Dataset %in% c("Test", "Validation")], aes(Prob, fill = Result)) +
  geom_histogram() +
  theme_bw(base_size = 16) +
  facet_grid(Dataset ~ Model, scales = "free_y") +
  theme(legend.position = "top") + 
  labs(x = "Probability", y = "Count")

ggsave("analysis/03.MachineLearning/Version7/02.Modeling/04.pcaNNet/02.CompareModels/Prediction_Probability_Distribution2.pdf", width = 16, height = 6)

