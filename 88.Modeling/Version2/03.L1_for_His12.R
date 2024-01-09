setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(caret)

outdir <- paste0("./analysis/81.ABFProcessing/SelectedSignals/01.StandardAA/06.Modeling/Version2")
if(!file.exists(outdir)) dir.create(outdir)

files <- list.files("./analysis/81.ABFProcessing/SelectedSignals/01.StandardAA/03.GroupedSignals", ".signal.txt", full.names = T)
Sigs <- lapply(files, fread)
Sigs <- data.table(A = rep(gsub(".signal.txt", "", basename(files)), mapply(nrow, Sigs)), do.call(rbind, Sigs))
Sigs[, ID0 := ID]
for (a in Sigs[, unique(AA)]) {
  Sigs[, ID0 := gsub(paste0(a, "_"), "", ID0)]
}
Sigs <- Sigs[AA %in% "His" & State == "Sington"]
FMfs <- list.files("./analysis/81.ABFProcessing/FeatureMatrix", "FeatureMatrix_", full.names = T, recursive = T)
FMfs <- unlist(lapply(Sigs[, unique(File)], FUN = function(x) grep(x, FMfs, value = T)))
FMs <- lapply(FMfs, fread)
FMs <- do.call(rbind, FMs)
if(!all(Sigs$ID0 %in% FMs$ID)) return(NULL)
FMs <- FMs[ID %in% Sigs$ID0]
setkey(FMs, ID)





set.seed(123)
trainIndex <- createDataPartition(Sigs[State == "Sington", AA], p = .8, 
                                  list = FALSE, 
                                  times = 1)
Sigs_Train <- Sigs[State == "Sington", ][trainIndex[, 1], ]
maxn <- pmin(Sigs_Train[, .N, A][, max(N)], 1000)

set.seed(234)
Sigs_Train_Upsample <- rbind(Sigs_Train[A %in% Sigs_Train[, .N, A][N < maxn, A], .SD[c(seq_len(.N), sample(.N, maxn - .N, replace = T))], A], 
                             Sigs_Train[A %in% Sigs_Train[, .N, A][N >= maxn, A], .SD[sample(.N, maxn, replace = F)], A])

setkey(Sigs_Train_Upsample, A, ID)

Sigs_Train_Upsample_FM <- as.data.frame(FMs[Sigs_Train_Upsample$ID0, grepl("X", colnames(FMs)) | colnames(FMs) %in% c("DeltaMean", "StageSD", "CurrentSD", "Segments", "SignalCurrentPercent", "DwellTime", "Blockade"), with = F])
Sigs_Train_Upsample_FM$Class <- factor(Sigs_Train_Upsample$A)

Sigs_Test <- Sigs[State == "Sington" & !ID %in% Sigs_Train_Upsample[, ID]]
Sigs_Test_FM <- as.data.frame(FMs[Sigs_Test$ID0, grepl("X", colnames(FMs)) | colnames(FMs) %in% c("DeltaMean", "StageSD", "CurrentSD", "Segments", "SignalCurrentPercent", "DwellTime", "Blockade"), with = F])

save(Sigs_Test, Sigs_Test_FM, Sigs_Train_Upsample, Sigs_Train_Upsample_FM, file = paste0(outdir, "/L1/His/Modeling_Data.RData"))



setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
outdir <- paste0("./analysis/81.ABFProcessing/SelectedSignals/01.StandardAA/06.Modeling/Version2")
load(paste0(outdir, "/L1/His/Modeling_Data.RData"))
library(data.table)
library(Biostrings)
library(caret)
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10)

library(doParallel)
cl <- makePSOCKcluster(20)
registerDoParallel(cl)

model1 <- train(Class ~ ., data = Sigs_Train_Upsample_FM, 
                # preProc = c("center", "scale", "YeoJohnson", "nzv"), 
                method = "rf", 
                trControl = fitControl,
                verbose = FALSE,
                ## to evaluate:
                # tuneGrid = expand.grid(mtry = 60),
                tuneLength = 6,
                metric = "Accuracy", 
                allowParallel = TRUE)
saveRDS(model1, file = paste0(outdir, "/L1/His/RF_model.Rds"))

stopCluster(cl)

