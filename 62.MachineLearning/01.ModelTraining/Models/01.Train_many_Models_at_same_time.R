setwd("/mnt/raid61/Personal_data/tangchao/AANanopore")
library(data.table)
library(Biostrings)
library(caret)
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10)

library(doParallel)
cl <- makePSOCKcluster(40)
registerDoParallel(cl)

ms <- c("LWQST", "LWQSTDE", "FYSL", "FYSLH", "FYSLHE", "LVFAG", "LVFA", "LVFAG", "LVFGKA", "LVFAK", "LVFGK", "LVFGK", "LVPFA", "LWQSTIFN", "LEF", "FHL", "LEFH", "LWQSTIFN", "LWQSTIFN", "LWQSTIFNDE", "LWQSTIFNDE", "LWQSTIFNDE", "DE", "DELWQSTIFN", "GAVLIFWYNQMSTP", "GAVLIFWYNQMSTPDE", "EKFS")
ms <- unique(ms)
ms <- strsplit(ms, "")
ms <- sort(mapply(function(x) paste(sort(x), collapse = ""), ms))
ms <- unique(ms)
ms <- ms[order(nchar(ms))]

lapply(ms, function(i) {
  AAs <- AMINO_ACID_CODE[strsplit(i, "")[[1]]]
  outdir <- paste0("./analysis/62.MachineLearning/01.ModelTraining/Model", paste0(sort(names(AMINO_ACID_CODE[AMINO_ACID_CODE %in% AAs])), collapse = ""))
  if(file.exists(paste0(outdir, "/model.Rds"))) {
    Fit1 <- readRDS(paste0(outdir, "/model.Rds"))
    if(identical(as.character(sort(Fit1$levels)), as.character(sort(AAs)))) return(NULL)
  }
  if(!dir.exists(outdir)) dir.create(outdir)
  
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
  if(!all(Sigs$ID0 %in% FMs$ID)) return(NULL)
  FMs <- FMs[ID %in% Sigs$ID0]
  setkey(FMs, ID)
  
  set.seed(123)
  trainIndex <- createDataPartition(Sigs$AA, p = .8, 
                                    list = FALSE, 
                                    times = 1)
  Sigs_Train <- Sigs[trainIndex[, 1], ]
  maxn <- Sigs_Train[, .N, AA][, max(N)]
  
  set.seed(234)
  Sigs_Train_Upsample <- rbind(Sigs_Train, Sigs_Train[, .SD[sample(.N, maxn - .N, replace = T)], AA])
  setkey(Sigs_Train_Upsample, AA, ID)
  
  Sigs_Train_Upsample_FM <- as.data.frame(FMs[Sigs_Train_Upsample$ID0, grep("X", colnames(FMs)), with = F])
  Sigs_Train_Upsample_FM$Class <- factor(Sigs_Train_Upsample$AA)
  
  Sigs_Test <- Sigs[!ID %in% Sigs_Train_Upsample[, ID]]
  Sigs_Test_FM <- as.data.frame(FMs[Sigs_Test$ID0, grep("X", colnames(FMs)), with = F])
  save(Sigs_Test, Sigs_Test_FM, Sigs_Train_Upsample, Sigs_Train_Upsample_FM, file = paste0(outdir,"/modeling_data.RData"))
  
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
})

stopCluster(cl)

