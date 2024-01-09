setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(caret)
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10)

library(doParallel)
cl <- makePSOCKcluster(20)
registerDoParallel(cl)



i <- paste(Biostrings::AA_STANDARD, collapse = "")

print(i)
AAs <- strsplit(i, "")[[1]]
AAs <- plyr::mapvalues(AAs, names(AMINO_ACID_CODE), AMINO_ACID_CODE, warn_missing = F)
if(is.element("Asx", AAs)) AAs <- gsub("Asx", "CbC", AAs)

outdir <- paste0("./analysis/81.ABFProcessing/SelectedSignals/01.StandardAA/03.TrainModels/Model", i, ".L2.Rds")
if(file.exists(outdir)) stop()

files <- list.files("./analysis/81.ABFProcessing/SelectedSignals/01.StandardAA/03.GroupedSignals", ".signal.txt", full.names = T)
Sigs <- do.call(rbind, lapply(files, fread))[State == "Mixtrue"]
Sigs[, ID0 := ID]
for (a in Sigs[, unique(AA)]) {
  Sigs[, ID0 := gsub(paste0(a, "_"), "", ID0)]
}
Sigs <- Sigs[AA %in% AAs]
FMfs <- list.files("./analysis/81.ABFProcessing/FeatureMatrix", "FeatureMatrix_", full.names = T, recursive = T)
FMfs <- unlist(lapply(Sigs[, unique(File)], FUN = function(x) grep(x, FMfs, value = T)))
FMs <- lapply(FMfs, fread)
FMs <- do.call(rbind, FMs)
if(!all(Sigs$ID0 %in% FMs$ID)) return(NULL)
FMs <- FMs[ID %in% Sigs$ID0]
setkey(FMs, ID)

set.seed(123)
trainIndex <- createDataPartition(Sigs$AA, p = .8, 
                                  list = FALSE, 
                                  times = 1)
Sigs_Train <- Sigs[trainIndex[, 1], ]
maxn <- pmin(Sigs_Train[, .N, AA][, max(N)], 3000)

set.seed(234)
Sigs_Train_Upsample <- rbind(Sigs_Train[AA %in% Sigs_Train[, .N, AA][N < maxn, AA], .SD[c(seq_len(.N), sample(.N, maxn - .N, replace = T))], AA], 
                             Sigs_Train[AA %in% Sigs_Train[, .N, AA][N >= maxn, AA], .SD[sample(.N, maxn, replace = F)], AA])

setkey(Sigs_Train_Upsample, AA, ID)

Sigs_Train_Upsample_FM <- as.data.frame(FMs[Sigs_Train_Upsample$ID0, grep("X", colnames(FMs)), with = F])
Sigs_Train_Upsample_FM$Class <- factor(Sigs_Train_Upsample$AA)

Sigs_Test <- Sigs[!ID %in% Sigs_Train_Upsample[, ID]]
Sigs_Test_FM <- as.data.frame(FMs[Sigs_Test$ID0, grep("X", colnames(FMs)), with = F])
save(Sigs_Test, Sigs_Test_FM, Sigs_Train_Upsample, Sigs_Train_Upsample_FM, file = gsub("Rds$", "RData", outdir))

set.seed(825)
Fit1 <- train(Class ~ ., data = Sigs_Train_Upsample_FM, 
              # preProc = c("center", "scale", "YeoJohnson", "nzv"),
              method = "rf", 
              trControl = fitControl,
              verbose = FALSE,
              ## to evaluate:
              tuneGrid = expand.grid(mtry = 2),
              # tuneLength = 2, 
              metric = "Accuracy", 
              allowParallel = TRUE)
saveRDS(Fit1, outdir)

stopCluster(cl)
