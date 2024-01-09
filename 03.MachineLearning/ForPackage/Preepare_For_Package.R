setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")

library(data.table)
library(ggplot2)
library(patchwork)
library(Biostrings)
library(cowplot)
library(parallel)
library(caret)
library(data.table)
library(Biostrings)
library(ggplot2)
library(parallel)

mean2 <- function(x, q1 = 0.25, q2 = 0.75) {
  q <- quantile(x, c(q1, q2))
  mean(x[x > min(q) & x < max(q)])
}
sd2 <- function(x, q1 = 0.1, q2 = 0.9) {
  q <- quantile(x, c(q1, q2))
  sd(x[x > min(q) & x < max(q)])
}

FeatureExtract <- function(x) {
  basemean <- x[L == "B", mean2(pA)]
  x[, pA := pA / basemean]
  dy <- round(density(x[L != "B", pA], from = 0, to = 1, n = 200, adjust = 1)$y, 3)
  names(dy) <- paste0("X", sprintf("%03d", 1:200))
  attr(dy, "AllTime") <- x[L != "B", diff(range(Time))] * 1000
  attr(dy, "DwellTime") <- x[Sm == x[L != "B", .N, Sm][which.max(N), Sm] & L != "B", diff(range(Time))] * 1000
  attr(dy, "SignalSD") <- x[L != "B", sd2(pA)]
  attr(dy, "Blockade") <- 1 - x[Sm == x[L != "B", .N, Sm][which.max(N), Sm] & L != "B", mean2(pA)]
  return(dy)
}

meta <- do.call(rbind, lapply(list.files("./analysis/01.AASignalRecognition/Version9/02.RangesL0L1/", ".txt", full.names = TRUE, recursive = T), fread))
meta <- meta[!file_name %in% c("21205012", "21205021", "21201008", "21303003")]
meta <- meta[, .(file_name, amino_acid, concentration, start_time, end_time)]
meta[, file_name := as.character(file_name)]

files <- list.files("./analysis/01.AASignalRecognition/Version9/04.FinalSignal", pattern = "_Sigs.txt", recursive = TRUE, full.names = T)
AAs <- lapply(files, fread)
names(AAs) <- gsub("_Sigs.txt", "", basename(files))
AAs <- AAs[!names(AAs) %in% c("21205012", "21205021", "21201008", "21303003")]

for(i in seq_along(AAs)) AAs[[i]] <- data.table(file_name = names(AAs)[i], ID = paste0(names(AAs)[i], "_", sprintf("%05d", seq_len(nrow(AAs[[i]])))), AAs[[i]])
AA <- do.call(rbind, AAs)

AA <- merge(AA, meta[, .(file_name, amino_acid)], by = "file_name")
setkey(AA, ID)



files <- list.files("./analysis/01.AASignalRecognition/Version9/03.BUBSignal", pattern = "_BUB.Rds", recursive = TRUE, full.names = T)
BUBs <- mclapply(files, readRDS, mc.cores = 10)
names(BUBs) <- gsub("_BUB.Rds", "", basename(files))
BUBs <- BUBs[!names(BUBs) %in% c("21205012", "21205021", "21201008", "21303003")]
for(i in seq_along(BUBs)) names(BUBs[[i]]) <- paste0(names(BUBs)[i], "_", sprintf("%05d", seq_along(BUBs[[i]])))
names(BUBs) <- NULL
BUBs <- do.call(c, BUBs)
BUBs <- BUBs[AA$ID]


Fs <- mclapply(BUBs, function(x) {
  FeatureExtract(x)
}, mc.cores = 10)

Fs_Tab <- mclapply(Fs, function(x) {
  t(rbind(data.frame(x = x, row.names = names(x)),
          data.frame(x = c(attr(x, "AllTime"), attr(x, "DwellTime"), attr(x, "SignalSD"), attr(x, "Blockade")),
                     row.names = c("AllTime", "DwellTime", "SignalSD", "Blockade"))))
}, mc.cores = 10)
Fs_Tab <- do.call(rbind, Fs_Tab)
row.names(Fs_Tab) <- AA$ID
Fs_Tab <- as.data.frame(Fs_Tab)

AA[, aa := plyr::mapvalues(amino_acid, AMINO_ACID_CODE, names(AMINO_ACID_CODE))]

Fs_Tab$Class <- factor(AA$aa)

Mat <- subset(Fs_Tab, AllTime > 4 & SignalSD < 4)
AAInfoi <- AA[ID %in% row.names(Mat)]


table(Mat$Class)

Mat <- data.table(Mat, keep.rownames = "ID")
set.seed(3456)
Train <- Mat[, .SD[sample(.N, round(Mat[, .N, Class][, min(N)]*0.8)), ], Class]
Test <- Mat[!ID %in% Train$ID]

Train <- data.frame(Train[, -c(1:2)], Class = Train[[1]])

save.image(file = "./analysis/03.MachineLearning/ForPackage/02.Modeling/01.RF/01.Models/RF_model.image.RData")
saveRDS(Train, file = "./analysis/03.MachineLearning/ForPackage/02.Modeling/01.RF/01.Models/Train.Rds")

fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10)

library(doParallel)
cl <- makePSOCKcluster(10)
registerDoParallel(cl)

RFC <- train(Class ~ ., 
               data = Train, 
               preProc = c("center", "scale", "YeoJohnson", "nzv"),
               method = "rf", 
               trControl = fitControl,
               verbose = FALSE,
               ## to evaluate:
               tuneGrid = expand.grid(mtry = 60),
               # tuneLength = 50,
               metric = "Accuracy", 
               allowParallel = TRUE)
saveRDS(RFC, file = paste0("./analysis/03.MachineLearning/ForPackage/02.Modeling/01.RF/01.Models/RF_model.Rds"))

saveRDS(RFC$finalModel, file = paste0("./analysis/03.MachineLearning/ForPackage/02.Modeling/01.RF/01.Models/RF_finalModel.Rds"))



set.seed(9560)
up_Mat <- upSample(x = Mat[, !colnames(Mat) %in% c("ID", "Class"), with = F], y = Mat$Class)
up_Mat <- as.data.table(up_Mat)
set.seed(3456)
Train <- up_Mat[, .SD[sample(.N, round(up_Mat[, .N, Class][, min(N)]*0.6)), ], Class]

saveRDS(Train, file = "./analysis/03.MachineLearning/ForPackage/02.Modeling/01.RF/01.Models/Train2.Rds")
Train <- readRDS("./analysis/03.MachineLearning/ForPackage/02.Modeling/01.RF/01.Models/Train2.Rds")




fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10)

library(doParallel)
cl <- makePSOCKcluster(40)
registerDoParallel(cl)

model <- train(Class ~ ., 
               data = Train, 
               preProc = c("center", "scale", "YeoJohnson", "nzv"),
               method = "rf", 
               trControl = fitControl,
               verbose = FALSE,
               ## to evaluate:
               tuneGrid = expand.grid(mtry = 60),
               # tuneLength = 50,
               metric = "Accuracy", 
               allowParallel = TRUE)
saveRDS(model, file = paste0("./analysis/03.MachineLearning/ForPackage/02.Modeling/01.RF/01.Models/RF_model_up.Rds"))
stopCluster(cl)

saveRDS(model$finalModel, file = paste0("./analysis/03.MachineLearning/ForPackage/02.Modeling/01.RF/01.Models/RF_model_up_finalModel.Rds"))



