setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(caret)

AAs <- c("Gln", "CbC", "Thr", "Ala", "Ile", "Phe", "Leu", "Ser", "Trp", "Arg", "His")
outdir <- paste0("./analysis/62.MachineLearning/01.ModelTraining/Model", paste0(names(AMINO_ACID_CODE[AMINO_ACID_CODE %in% AAs]), collapse = ""))
dir.create(outdir)

files <- list.files("./analysis/61.SignalSelecting/01.StandardAA/02.SelectedSignals", "_V2.Rds", full.names = T)
Sigs <- lapply(files, function(x) {
  readRDS(x)[[1]]
})
names(Sigs) <- gsub("_State_Signals_V2.Rds", "", basename(files))
Sigs <- data.table(AA = rep(names(Sigs), mapply(nrow, Sigs)), do.call(rbind, Sigs))
Sigs <- Sigs[State == "Sington"]
Sigs <- Sigs[A %in% AAs]
Sigs[, ID0 := ID]
for (a in Sigs[, unique(A)]) {
  Sigs[, ID0 := gsub(paste0(a, "_"), "", ID0)]
}

FMfs <- list.files("./analysis/21.ABFProcessing", "FeatureMatrix_", full.names = T, recursive = T)
FMfs <- mapply(Sigs[, unique(File)], FUN = function(x) grep(x, FMfs, value = T)[1])
FMs <- lapply(FMfs, fread)
FMs <- do.call(rbind, FMs)
FMs <- FMs[ID %in% Sigs$ID0]
setkey(FMs, ID)

set.seed(123)
trainIndex <- createDataPartition(Sigs$A, p = .8, 
                                  list = FALSE, 
                                  times = 1)
Sigs_Train <- Sigs[trainIndex[, 1], ]
maxn <- Sigs_Train[, .N, A][, max(N)]

set.seed(234)
Sigs_Train_Upsample <- rbind(Sigs_Train, Sigs_Train[, .SD[sample(.N, maxn - .N, replace = T)], A])
setkey(Sigs_Train_Upsample, A, ID)

Sigs_Train_Upsample_FM <- as.data.frame(FMs[Sigs_Train_Upsample$ID0, grep("X", colnames(FMs)), with = F])
Sigs_Train_Upsample_FM$Class <- factor(Sigs_Train_Upsample$A)

Sigs_Test <- Sigs[!ID %in% Sigs_Train_Upsample[, ID]]
Sigs_Test_FM <- as.data.frame(FMs[Sigs_Test$ID0, grep("X", colnames(FMs)), with = F])

save(Sigs_Test, Sigs_Test_FM, Sigs_Train_Upsample, Sigs_Train_Upsample_FM, file = paste0(outdir,"/modeling_data.RData"))

fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10)

library(doParallel)
cl <- makePSOCKcluster(20)
registerDoParallel(cl)

# library(doSNOW)
# cl <- makeSOCKcluster(2)
# registerDoSNOW(cl)

set.seed(825)
Fit1 <- train(Class ~ ., data = Sigs_Train_Upsample_FM, 
              preProc = c("center", "scale", "YeoJohnson", "nzv"),
              method = "rf", 
              trControl = fitControl,
              verbose = FALSE,
              ## to evaluate:
              tuneGrid = expand.grid(mtry = 2),
              # tuneLength = 2, 
              metric = "Accuracy", 
              allowParallel = TRUE)
saveRDS(Fit1, paste0(outdir,"/model.Rds"))
stopCluster(cl)




