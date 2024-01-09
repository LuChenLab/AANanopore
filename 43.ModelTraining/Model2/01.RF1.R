setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(shiny)
library(plotly)
library(Biostrings)

set.seed(19)
AA_Cols <- sample(c(RColorBrewer::brewer.pal(n = 12, "Paired"), RColorBrewer::brewer.pal(n = 8, "Dark2")), 20)
names(AA_Cols) <- Biostrings::AA_STANDARD
AA_Cols[AA_Cols == "#FFFF99"] <- "#1A1A1A"

get_density <- function(x, y, n = 1000, ...) {
  dens <- MASS::kde2d(x, y, n = n, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii] / max(dens$z[ii]))
}


AABlockade <- data.table(AA = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"), 
                         Blockade = c(0.14718, 0.17148, 0.16516, 0.21409, 0.16, 0.24528, 0.1882, 0.12072, 0.24652, 0.20722, 0.1995, 0.16875, 0.19772, 0.22018, 0.2183, 0.13131, 0.16101, 0.22744, 0.21276, 0.19044))
AABlockade$amino_acid <- plyr::mapvalues(AABlockade$AA, names(Biostrings::AMINO_ACID_CODE), Biostrings::AMINO_ACID_CODE)

AAGroup <- list(GSA = c('G', 'S', 'A'), 
                TNRK = c("T", "N", "R", "K"), 
                QVML = c("Q", "V", "M", "L"), 
                IYDPFW = c("I", "Y", "D", "P", "F", "W"), 
                EH = c("E", "H"))
AAGroup <- data.table(Group = rep(names(AAGroup), mapply(length, AAGroup)), AA = unlist(AAGroup))
AAGroup[, amino_acid := plyr::mapvalues(AA, names(AMINO_ACID_CODE), AMINO_ACID_CODE)]


Sigs <- lapply(AAGroup[, amino_acid], function(aat) {
  sigs <- do.call(rbind, lapply(list.files(paste0("./analysis/42.SignalSelecting/01.StandardAA/", aat), full.names = T), fread))
  sigs$File <- stringr::str_remove_all(sigs$ID, "_([[:digit:]]+)$")
  sigs$D <- sigs[, get_density(x = Blockade, y = log10(DwellTime))]
  sigs
})
Sigs <- do.call(rbind, Sigs)
Sigs <- merge(Sigs, AAGroup, by.x = "A", by.y = "amino_acid")

ggplot(Sigs[Group == "GSA" & D > .1], aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_point() + 
  scale_y_log10()
ggplot(Sigs[Group == "GSA" & D > .1], aes(x = Blockade, fill = AA)) + 
  geom_histogram(binwidth = 0.0003)


ggplot(Sigs[Group == "TNRK" & D > .1], aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_point() + 
  scale_y_log10()
ggplot(Sigs[Group == "TNRK" & D > .1], aes(x = Blockade, fill = AA)) + 
  geom_histogram(binwidth = 0.0003)


ggplot(Sigs[Group == "QVML" & D > .1], aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_point() + 
  scale_y_log10()
ggplot(Sigs[Group == "QVML" & D > .1], aes(x = Blockade, fill = AA)) + 
  geom_histogram(binwidth = 0.0003)



ggplot(Sigs[Group == "IYDPFW" & D > .1], aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_point() + 
  scale_y_log10()
ggplot(Sigs[Group == "IYDPFW" & D > .1], aes(x = Blockade, fill = AA)) + 
  geom_histogram(binwidth = 0.0003)


ggplot(Sigs[Group == "EH" & D > .1], aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_point() + 
  scale_y_log10()
ggplot(Sigs[Group == "EH" & D > .1], aes(x = Blockade, fill = AA)) + 
  geom_histogram(binwidth = 0.0003)


ggplot(Sigs[D > .15], aes(x = Blockade, y = DwellTime, colour = AA)) + 
  geom_point(alpha = .5) + 
  scale_y_log10() + 
  facet_wrap(~ Group, scales = 'free_x') + 
  scale_colour_manual(breaks = names(AA_Cols), values = AA_Cols)

ggplot(Sigs[D > .15], aes(x = Blockade, fill = AA, colour = AA)) + 
  geom_density(alpha = .5, adjust = 2) + 
  facet_wrap(~ Group, scales = 'free_x') + 
  scale_fill_manual(breaks = names(AA_Cols), values = AA_Cols) + 
  scale_colour_manual(breaks = names(AA_Cols), values = AA_Cols)


Sigs[D > .15, .N, AA][order(N)]
Sigs[D > .1, .N, AA][order(N)]


SigsID <- Sigs[D > .1, .(ID, AA)]

NoisID <- fread("analysis/42.SignalSelecting/02.Background/Noise_Signal_ID.txt")

IDTU <- rbind(SigsID, NoisID)
IDTU[, AA := factor(AA, levels = setdiff(c(AA_STANDARD, 'Noise'), 'C'))]
IDTU[, .N, AA]

set.seed(123)
Test_Set <- IDTU[, .SD[sample(.N, floor(.N * .2))], AA]
Train_Set <- IDTU[!ID %in% Test_Set[, ID]]
Train_Set[, .N, AA][order(N)]


Selected_Sig_Upsample <- rbind(Train_Set[AA %in% c('N', 'Noise'), .SD[sample(.N, 3500)], AA], 
                               Train_Set[!AA %in% c('N', 'Noise'), .SD[sample(.N, 3500 - .N, replace = T)], AA], 
                               Train_Set[!AA %in% c('N', 'Noise'), ])
Selected_Sig_Upsample[, AA := factor(AA, levels = setdiff(c(AA_STANDARD, 'Noise'), 'C'))]
setkey(Selected_Sig_Upsample, AA, ID)

files <- list.files("./analysis/21.ABFProcessing/02.Background/FeatureMatrix", pattern = "", full.names = TRUE, recursive = TRUE)
Feature_Mat0 <- lapply(files, fread)
Feature_Mat0 <- do.call(rbind, Feature_Mat0)
Feature_Mat0 <- na.omit(Feature_Mat0)

files <- list.files("./analysis/21.ABFProcessing/01.StandardAA/FeatureMatrix", pattern = "", full.names = TRUE, recursive = TRUE)
Feature_Mat <- lapply(files, fread)
Feature_Mat <- do.call(rbind, Feature_Mat)

setkey(Feature_Mat0, ID)
Feature_Mat <- rbind(Feature_Mat, Feature_Mat0)
setkey(Feature_Mat, ID)



Train <- setDF(Feature_Mat[Selected_Sig_Upsample[, ID], ])
Train <- Train[, c(paste0("X", sprintf("%04d", 601:1000)), paste0("P", sprintf("%03d", 1:100)), "StageSD", "CurrentSD", "Segments", "SignalCurrentPercent", "SignalCurrentWidth", "Blockade", "DwellTime")]
Train$Class <- Selected_Sig_Upsample[, AA]

Test <- setDF(Feature_Mat[Test_Set[, ID], ])
Test <- Test[, c(paste0("X", sprintf("%04d", 601:1000)), paste0("P", sprintf("%03d", 1:100)), "StageSD", "CurrentSD", "Segments", "SignalCurrentPercent", "SignalCurrentWidth", "Blockade", "DwellTime")]
Test$Class <- Test_Set[, AA]


saveRDS(Train, "./analysis/43.ModelTraining/Model2/Train1.Rds")
saveRDS(Test, "./analysis/43.ModelTraining/Model2/Test1.Rds")


setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(ggplot2)
library(parallel)
library(caret)

fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10)

library(doParallel)
cl <- makePSOCKcluster(20)
registerDoParallel(cl)

Train <- readRDS("./analysis/43.ModelTraining/Model2/Train1.Rds")
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
saveRDS(RFC, file = "./analysis/43.ModelTraining/Model2/RF1.Rds")

stopCluster(cl)



RFC <- readRDS("./analysis/43.ModelTraining/Model2/RF1.Rds")

Res <- data.table(Obse = Test$Class, Pred = predict(RFC, Test), Prob = apply(predict(RFC, Test, type = "prob"), 1, max))
hist(Res[, Prob], breaks = 100)

Res[, Obse := factor(Obse, levels = c('G', 'S', 'A', 'T', 'N', 'R', 'K', 'V', 'Q', 'L', 'M', 'I', 'Y', 'D', 'P', 'F', 'W', 'E', 'H', 'Noise'))]
Res[, Pred := factor(Pred, levels = c('G', 'S', 'A', 'T', 'N', 'R', 'K', 'V', 'Q', 'L', 'M', 'I', 'Y', 'D', 'P', 'F', 'W', 'E', 'H', 'Noise'))]
cfM <- as.matrix(Res[, confusionMatrix(Obse, Pred)])
cfM <- t(t(cfM)/colSums(cfM))
pheatmap::pheatmap(cfM, cluster_rows = F, cluster_cols = F, display_numbers = T)


Train_Set[AA %in% c('G', 'Noise'), .N, AA]
Train_Set[AA %in% c('T', 'N', 'R', 'K'), .N, AA]
Train_Set[AA %in% c('Q', 'V'), .N, AA]
Train_Set[AA %in% c('L', 'M'), .N, AA]
Train_Set[AA %in% c('I', 'Y', 'D', 'P', 'F', 'W'), .N, AA]

TrainSub1 <- subset(Train, Class %in% c('G', 'Noise'))
TrainSub2 <- subset(Train, Class %in% c('T', 'N', 'R', 'K'))
TrainSub3 <- subset(Train, Class %in% c('Q', 'V'))
TrainSub4 <- subset(Train, Class %in% c('L', 'M'))
TrainSub5 <- subset(Train, Class %in% c('I', 'Y', 'D', 'P', 'F', 'W'))

TrainSub1$Class <- factor(TrainSub1$Class, levels = c('G', 'Noise'))
TrainSub2$Class <- factor(TrainSub2$Class, levels = c('T', 'N', 'R', 'K'))
TrainSub3$Class <- factor(TrainSub3$Class, levels = c('Q', 'V'))
TrainSub4$Class <- factor(TrainSub4$Class, levels = c('L', 'M'))
TrainSub5$Class <- factor(TrainSub5$Class, levels = c('I', 'Y', 'D', 'P', 'F', 'W'))


Train <- readRDS("./analysis/43.ModelTraining/Model2/Train1.Rds")
Train <- subset(Train, Class %in% c('G', 'Noise'))
Train$Class <- factor(Train$Class, levels = c('G', 'Noise'))
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
saveRDS(RFC, file = "./analysis/43.ModelTraining/Model2/RF1_1.Rds")




Train <- readRDS("./analysis/43.ModelTraining/Model2/Train1.Rds")
Train <- subset(Train, Class %in% c('T', 'N', 'R', 'K'))
Train$Class <- factor(Train$Class, levels = c('T', 'N', 'R', 'K'))
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
saveRDS(RFC, file = "./analysis/43.ModelTraining/Model2/RF1_2.Rds")




Train <- readRDS("./analysis/43.ModelTraining/Model2/Train1.Rds")
Train <- subset(Train, Class %in% c('Q', 'V'))
Train$Class <- factor(Train$Class, levels = c('Q', 'V'))
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
saveRDS(RFC, file = "./analysis/43.ModelTraining/Model2/RF1_3.Rds")




Train <- readRDS("./analysis/43.ModelTraining/Model2/Train1.Rds")
Train <- subset(Train, Class %in% c('L', 'M'))
Train$Class <- factor(Train$Class, levels = c('L', 'M'))
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
saveRDS(RFC, file = "./analysis/43.ModelTraining/Model2/RF1_4.Rds")




Train <- readRDS("./analysis/43.ModelTraining/Model2/Train1.Rds")
Train <- subset(Train, Class %in% c('I', 'Y', 'D', 'P', 'F', 'W'))
Train$Class <- factor(Train$Class, levels = c('I', 'Y', 'D', 'P', 'F', 'W'))
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
saveRDS(RFC, file = "./analysis/43.ModelTraining/Model2/RF1_5.Rds")
